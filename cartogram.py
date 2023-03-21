import pyfftw

import border_util

import geopandas as gpd
import numpy as np

from shapely.affinity import scale, translate
from shapely.geometry import box


def paired_distance(X, Y):
    Z = X - Y
    norms = np.einsum("ij, ij -> i", Z, Z)
    return np.sqrt(norms, norms)


class Cartogram:
    def __init__(self, gdf, value_field, id_field=None, geometry_field='geometry'):

        # Initialise Gastner-Newman variables
        self.xsize = None
        self.ysize = None
        self.rhot = None
        self.fftrho = None
        self.fftexpt = None
        self.expky = None
        self.grid_points_x = None
        self.grid_points_y = None
        self.errorp = None
        self.drp = None
        self.inith = 0.001
        self.target_error = 0.01
        self.max_ratio = 4.0
        self.expected_time = 100000000.0

        self.gdf = gdf
        self.value_field = value_field
        self.geo_field = geometry_field

        if not id_field:
            self.gdf['id_field'] = self.gdf.index
            self.id_field = "id_field"
        else:
            self.id_field = id_field


    def non_contiguous(self, position="centroid", size_value=1.0):
        geodf = self.gdf[[self.value_field, self.id_field, self.geo_field]].copy()

        geo = geodf[self.geo_field]
        positions = {
            "centroid": geo.centroid,
            "centre": geo.envelope.centroid
        }

        if position.lower() in positions.keys():
            geodf["cent"] = positions[position.lower()]
        else:
            print("Incorrect Position Parameter, use: 'centroid' | 'centre'")

        geodf["density"] = geodf[self.value_field] / geodf.area
        geodf["rank"] = geodf["density"].rank(axis=0, method="first", ascending=False)

        anchor = geodf[geodf["rank"] == 1]["density"].values[0]

        geodf["scale"] = (1.0 / np.power(anchor, 0.5)) * np.power(geodf[self.value_field] / geodf.area, 0.5) * size_value

        new_geo = [
            scale(
                g[1][self.geo_field],
                xfact=g[1]["scale"],
                yfact=g[1]["scale"],
                origin=g[1]["cent"]
            )
            for g in geodf.iterrows()
        ]

        del geodf["density"], geodf["rank"], geodf["cent"]

        return gpd.GeoDataFrame(geodf, geometry=new_geo)


    def dorling(self, ratio=0.4, friction=0.25, iterations=100, stop=None):
        def repel(neighbour, focal, xrepel, yrepel):
            xrepel -= (
                neighbour["overlap"] * (neighbour["geometry"].x - focal["geometry"].x) / neighbour["dist"]
            )
            yrepel -= (
                neighbour["overlap"] * (neighbour["geometry"].y - focal["geometry"].y) / neighbour["dist"]
            )

            return xrepel, yrepel

        def attract(nb, borders, idx, focal, perimeter, xattract, yattract):
            if sum((borders["focal"] == idx) & (borders["neighbor"] == nb.name)) == 1:
                nb["overlap"] = (
                    np.abs(nb["overlap"])
                    * float(
                        borders[(borders["focal"] == idx) & (borders["neighbor"] == nb.name)]["weight"]
                    )
                    / perimeter[idx]
                )

            xattract += nb["overlap"] * (nb["geometry"].x - focal["geometry"].x) / nb["dist"]
            yattract += nb["overlap"] * (nb["geometry"].y - focal["geometry"].y) / nb["dist"]

            return xattract, yattract

        borders, islands = border_util.get_borders(self.gdf)
        perimeter = self.gdf.length

        regions = gpd.GeoDataFrame(
            self.gdf.drop(columns=self.geo_field), geometry=self.gdf.centroid
        )

        focal = np.stack(
            borders.merge(
                regions[self.geo_field].map(np.array).to_frame(),
                left_on="focal",
                right_index=True,
            ).sort_index()[self.geo_field]
        )

        neighbour = np.stack(
            borders.merge(
                regions[self.geo_field].map(np.array).to_frame(),
                left_on="neighbor",
                right_index=True,
            ).sort_index()[self.geo_field]
        )

        total_distance = np.sum(paired_distance(focal, neighbour))

        focal_radius = borders.merge(
            regions[[self.value_field]],
            left_on="focal",
            right_index=True,
        ).sort_index()[self.value_field]

        neighbour_radius = borders.merge(
            regions[[self.value_field]],
            left_on="neighbor",
            right_index=True,
        ).sort_index()[self.value_field]

        total_radius = np.sum(
            (focal_radius / np.pi) ** 0.5 + (neighbour_radius / np.pi) ** 0.5
        )

        scale = total_distance / total_radius

        regions["radius"] = np.power(regions[self.value_field] / np.pi, 0.5) * scale
        widest = regions["radius"].max()

        for i in range(iterations):
            print(f"Starting Iteration: {i}")
            displacement = 0.0

            for idx, region in regions.iterrows():
                if stop is not None:
                    if idx == stop + 1:
                        break
                xrepel = 0.0
                yrepel = 0.0
                xattract = 0.0
                yattract = 0.0
                closest = widest

                neighbours = regions[
                    regions.distance(region[self.geo_field]).between(
                        0, widest + region["radius"], inclusive="neither",
                    )
                ].copy()

                if len(neighbours) > 0:
                    neighbours["dist"] = neighbours[self.geo_field].distance(region[self.geo_field])

                    closest = widest if neighbours["dist"].min() > widest else neighbours["dist"].min()

                    neighbours["overlap"] = (neighbours["radius"] + region["radius"]) - neighbours["dist"]

                    for idy, nb in neighbours.iterrows():
                        if nb["overlap"] > 0.0:
                            xrepel, yrepel = repel(nb, region, xrepel, yrepel)
                        else:
                            xattract, yattract = attract(nb, borders, idx, region, perimeter, xattract, yattract)

                attract_dist = np.sqrt((xattract ** 2) + (yattract ** 2))
                repel_dist = np.sqrt((xrepel ** 2) + (yrepel ** 2))

                if repel_dist > closest:
                    xrepel = closest * xrepel / (repel_dist + 1.0)
                    yrepel = closest * yrepel / (repel_dist + 1.0)
                    repel_dist = closest

                if repel_dist > 0:
                    xtotal = (1.0 - ratio) * xrepel + ratio * (
                        repel_dist * xattract / (attract_dist + 1.0)
                    )
                    ytotal = (1.0 - ratio) * yrepel + ratio * (
                        repel_dist * yattract / (attract_dist + 1.0)
                    )
                else:
                    if attract_dist > closest:
                        xattract = closest * xattract / (attract_dist + 1.0)
                        yattract = closest * yattract / (attract_dist + 1.0)
                    xtotal = xattract
                    ytotal = yattract

                displacement += np.sqrt((xtotal ** 2) + (ytotal ** 2))

                xvector = friction * xtotal
                yvector = friction * ytotal

                regions.loc[idx, self.geo_field] = translate(
                    region[self.geo_field], xoff=xvector, yoff=yvector
                )

            displacement = displacement / len(regions)

        return gpd.GeoDataFrame(
            data=regions.drop(columns=["geometry", "radius"]),
            geometry=regions.apply(lambda x: x["geometry"].buffer(x["radius"]), axis=1)
        )

    def diffusion(self, max_iter=100, diff_coeff=0.25, grid_size=(100, 100)):
        def make_density_grid():
            minx, miny, maxx, maxy = self.gdf.total_bounds

            gdf = self.gdf.copy()
            gdf['density'] = gdf[self.value_field] / gdf.area
            mean_density = gdf['density'].mean()

            W = maxx - minx
            H = maxy - miny

            n_cells_x, n_cells_y = grid_size

            cell_x = W / n_cells_x
            cell_y = H / n_cells_y

            padding = n_cells_x

            n_cells_x_padded = n_cells_x + 2 * padding
            n_cells_y_padded = n_cells_y + 2 * padding

            minx_padded = minx - (padding * cell_x)
            miny_padded = miny - (padding * cell_y)
            maxx_padded = maxx + (padding * cell_x)
            maxy_padded = maxy + (padding * cell_y)

            x_coords_padded = np.arange(minx_padded, maxx_padded, cell_x)
            y_coords_padded = np.arange(miny_padded, maxy_padded - cell_y, cell_y)

            geometries = [box(x, y, x + cell_x, y + cell_y) for y in y_coords_padded for x in x_coords_padded]
            fishnet = gpd.GeoDataFrame(geometry=geometries, crs=gdf.crs)

            sindex = gdf.sindex

            densities = []

            for idx, cell in fishnet.iterrows():
                cell_geom = cell['geometry']

                possible_matches_idx = list(sindex.intersection(cell_geom.bounds))
                possible_matches = gdf.iloc[possible_matches_idx]

                precise_matches = possible_matches[possible_matches.intersects(cell_geom)]

                if not precise_matches.empty:
                    density = precise_matches['density'].max()
                else:
                    density = mean_density

                densities.append(density)

            fishnet['density'] = densities

            return fishnet, n_cells_x_padded, n_cells_y_padded


        def create_grid_of_points():
            self.grid_points_x = np.zeros((self.xsize + 1) * (self.ysize + 1))
            self.grid_points_y = np.zeros((self.xsize + 1) * (self.ysize + 1))

            i = 0

            for iy in range(0, self.ysize+1):
                for ix in range(0, self.xsize+1):
                    self.grid_points_x[i] = ix
                    self.grid_points_y[i] = iy
                    i += 1

        def compute_initial_density(density_grid):
            densities = density_grid['density'].values

            # Convert density grid to column-major form
            rho = np.array(densities).reshape((self.xsize, self.ysize), order='F')

            self.fftrho = pyfftw.empty_aligned((self.xsize, self.ysize), dtype='complex128')

            # TODO: rfftn and fftn do not produce expected outputs; try changing to scipy's DCT

            fft_plan = pyfftw.builders.rfftn(self.fftrho, s=(self.xsize, self.ysize), axes=(0, 1), planner_effort='FFTW_ESTIMATE', norm='forward')
            self.fftrho = fft_plan(rho)


        def density_snapshot(t, s):
            kx, ky = 0, 0
            expkx = 0

            # Calculate expky array to save time in next part
            for iy in range(0, self.ysize):
                ky = np.pi * iy / self.ysize
                self.expky[iy] = np.exp(-ky * ky * t)

            for ix in range(0, self.xsize):
                kx = np.pi * ix / self.xsize
                expkx = np.exp(-kx * kx * t)

                for iy in range (0, self.ysize):
                    self.fftexpt[ix][iy] = expkx * self.expky[iy] * self.fftrho[ix][iy]
                    self.rhot[s][ix][iy] = self.fftexpt[ix][iy]

            # Perform back-transform
            fft_plan = pyfftw.builders.irfftn(self.rhot[s], s=(self.xsize, self.ysize), axes=(0, 1), planner_effort='FFTW_ESTIMATE')
            self.fftrho = fft_plan(self.rhot[s])


        def velocity(rx, ry, s):
            ix = int(rx)

            if ix < 0:
                ix = 0
            elif ix >= self.xsize:
                ix = self.xsize - 1

            ixm1 = ix - 1
            if ixm1 < 0:
                ixm1 = 0

            ixp1 = ix + 1
            if ixp1 >= self.xsize:
                ixp1 = self.xsize - 1

            iy = int(ry)
            if iy < 0:
                iy = 0
            elif iy >= self.ysize:
                iy = self.ysize - 1

            iym1 = iy - 1
            if iym1 < 0:
                iym1 = 0

            iyp1 = iy + 1
            if iyp1 >= self.ysize:
                iyp1 = self.ysize - 1

            # Calculate densities at nine surrounding grid points
            rho00 = self.rhot[s][ixm1][iym1]
            rho10 = self.rhot[s][ix][iym1]
            rho20 = self.rhot[s][ixp1][iym1]
            rho01 = self.rhot[s][ixm1][iy]
            rho11 = self.rhot[s][ix][iy]
            rho21 = self.rhot[s][ixp1][iy]
            rho02 = self.rhot[s][ixm1][iyp1]
            rho12 = self.rhot[s][ix][iyp1]
            rho22 = self.rhot[s][ixp1][iyp1]


            # Calculate velocities at four surrounding grid points
            mid11 = rho00 + rho10 + rho01 + rho11
            vx11 = -2.0 * (rho10 - rho00 + rho11 - rho01) / mid11
            vy11 = -2.0 * (rho01 - rho00 + rho11 - rho10) / mid11

            mid21 = rho10 + rho20 + rho11 + rho21
            vx21 = -2.0 * (rho20 - rho10 + rho21 - rho11) / mid21
            vy21 = -2.0 * (rho11 - rho10 + rho21 - rho20) / mid21

            mid12 = rho01 + rho11 + rho02 + rho12
            vx12 = -2.0 * (rho11 - rho01 + rho12 - rho02) / mid12
            vy12 = -2.0 * (rho02 - rho01 + rho12 - rho11) / mid12

            mid22 = rho11 + rho21 + rho12 + rho22
            vx22 = -2.0 * (rho21 - rho11 + rho22 - rho12) / mid22
            vy22 = -2.0 * (rho12 - rho11 + rho22 - rho21) / mid22


            # Calculate weights for bilinear interpolation
            dx = rx - ix
            dy = ry - iy

            dx1m = 1.0 - dx
            dy1m = 1.0 - dy

            w11 = dx1m * dy1m
            w21 = dx * dy1m
            w12 = dx1m * dy
            w22 = dx * dy

            # Perform interpolation for x and y components of velocity
            vxp = w11*vx11 + w21*vx21 + w12*vx12 + w22*vx22
            vyp = w11*vy11 + w21*vy21 + w12*vy12 + w22*vy22

            vp = (vxp, vyp)

            return vp


        def integrate_two_steps(t, h, s):
            s0 = s
            s1 = (s + 1) % 5
            s2 = (s + 2) % 5
            s3 = (s + 3) % 5
            s4 = (s + 4) % 5

            density_snapshot((t + 0.5 * h), s1)
            density_snapshot((t + 1.0 * h), s2)
            density_snapshot((t + 1.5 * h), s3)
            density_snapshot((t + 2.0 * h), s4)


            # Do three Runga-Kutta steps for each point in turn
            esq_max = 0.0
            drsq_max = 0.0
            n_points = (self.xsize + 1) * (self.ysize + 1)

            for p in range(0, n_points):
                rx1 = self.grid_points_x[p]
                ry1 = self.grid_points_y[p]

                # Do big combined 2-Hour Runga-Kutta step
                v1 = velocity(rx1, ry1, s0)
                k1x = 2 * h * v1[0]
                k1y = 2 * h * v1[1]

                v2 = velocity((rx1 + 0.5*k1x), (ry1 + 0.5*k1y), s2)
                k2x = 2 * h * v2[0]
                k2y = 2 * h * v2[1]

                v3 = velocity((rx1 + 0.5*k2x), (ry1 + 0.5*k1y), s3)
                k3x = 2 * h * v3[0]
                k3y = 2 * h * v3[1]

                v4 = velocity((rx1 + 0.5*k3x), (ry1 + 0.5*k3y), s3)
                k4x = 2 * h * v4[0]
                k4y = 2 * h * v4[1]

                dx12 = (k1x + k4x + 2.0 * (k2x + k3x)) / 6.0
                dy12 = (k1y + k4y + 2.0 * (k2y + k3y)) / 6.0

                # Do first small RK step
                k1x = h * v1[0]
                k1y = h * v1[1]

                v2 = velocity((rx1 + 0.5 * k1x), (ry1 + 0.5 * k1y), s1)
                k2x = h * v2[0]
                k2y = h * v2[1]

                v3 = velocity((rx1 + 0.5 * k2x), (ry1 + 0.5 * k2y), s1)
                k3x = h * v3[0]
                k3y = h * v3[1]

                v4 = velocity((rx1 + k3x), (ry1 + k3y), s2)
                k4x = h * v4[0]
                k4y = h * v4[1]

                dx1 = (k1x + k4x + 2.0 * (k2x + k3x)) / 6.0
                dy1 = (k1y + k4y + 2.0 * (k2y + k3y)) / 6.0

                # Do second small RK step
                rx2 = rx1 + dx1
                ry2 = ry1 + dy1

                v1 = velocity(rx2, ry2, s2)
                k1x = h * v1[0]
                k1y = h * v1[1]

                v2 = velocity((rx2 + 0.05*k1x), (ry2 + 0.5*k1y), s3)
                k2x = h * v2[0]
                k2y = h * v2[1]

                v3 = velocity((rx2 + 0.5*k2x), (ry2 + 0.5*k2y), s3)
                k3x = h * v3[0]
                k3y = h * v3[1]

                v4 = velocity((rx2 + k3x), (ry2 + k3y), s4)
                k4x = h * v4[0]
                k4y = h * v4[1]

                dx2 = (k1x + k4x + 2.0 * (k2x + k3x)) / 6.0
                dy2 = (k1y + k4y + 2.0 * (k2y + k3y)) / 6.0

                # Calculate Squared Error
                ex = (dx1 + dx2 - dx12) / 15
                ey = (dy1 + dy2 - dy12) / 15
                esq = ex*ex + ey*ey

                if esq > esq_max:
                    esq_max = esq

                # Update position of vertex using more accurate result
                # and deal with Boundary Conditions
                dxtotal = dx1 + dx2 + ex
                dytotal = dy1 + dy2 + ey

                drsq = dxtotal*dxtotal + dytotal*dytotal

                if drsq > drsq_max:
                    drsq_max = drsq

                rx3 = rx1 + dxtotal
                ry3 = ry1 + dytotal

                if rx3 < 0:
                    rx3 = 0
                elif rx3 > self.xsize:
                    rx3 = self.xsize

                if ry3 < 0:
                    ry3 = 0
                elif ry3 > self.ysize:
                    ry3 = self.ysize

                self.grid_points_x[p] = rx3
                self.grid_points_y[p] = ry3

            self.errorp = np.sqrt(esq_max)
            self.drp = np.sqrt(drsq_max)

            return s4


        def update_running_status(t):
            perc = int(np.round(100.0 * np.log(t / self.inith) / np.log(self.expected_time / self.inith)))

            if perc > 100:
                perc = 100

            res = (700 - 350) * perc / 100
            res += 350

            print(f"Diffusion Process: {perc}% done")

        def make_cartogram(blur):
            s = 0
            density_snapshot(0.0, s)

            step = 0
            t = 0.5 * blur * blur
            h = self.inith

            sp = None
            self.drp = 1.0

            while self.drp > 0.0:
                sp = integrate_two_steps(t, h, s)

                self.drp = 0.0

                t += 2.0 * h
                step += 2
                s = sp

                # Adjust time-step
                desired_ratio = np.power((2 * self.target_error / self.errorp), 0.2)

                if desired_ratio > self.max_ratio:
                    h *= self.max_ratio
                else:
                    h *= desired_ratio

                update_running_status(t)


        def project_cartogram():
            pass

        density_grid, self.xsize, self.ysize = make_density_grid()

        densities = density_grid['density'].values

        # Convert density grid to column-major form
        # rho = np.array(densities).reshape((self.xsize, self.ysize), order='F')

        self.rhot = np.zeros((5, self.xsize, self.ysize))
        self.fftrho = np.zeros((self.xsize, self.ysize))
        self.fftexpt = np.zeros((self.xsize, self.ysize))
        self.expky = np.zeros(self.ysize)

        compute_initial_density(density_grid)

        create_grid_of_points()

        make_cartogram(0.0)

        # transform = project_cartogram(cartogram)

        return density_grid

    def fast_flow(self):
        pass