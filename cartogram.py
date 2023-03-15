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

            x_coords = np.arange(minx, maxx, cell_x)
            y_coords = np.arange(miny, maxy, cell_y)

            geometries = [box(x, y, x + cell_x, y + cell_y) for y in y_coords for x in x_coords]
            fishnet = gpd.GeoDataFrame(geometry=geometries, crs=gdf.crs)

            sindex = gdf.sindex

            densities = []

            for i, cell in fishnet.iterrows():
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

            return fishnet

    # TODO: Implement rest of Gastner-Newman diffusion algorithm

    def fast_flow(self):
        pass