from functools import partial

import pyfftw
import geopandas as gpd
import numpy as np
import pandas as pd

from shapely.affinity import scale, translate
from libpysal.weights import Rook, W, KNN, attach_islands


def paired_distance(X, Y):
    Z = X - Y
    norms = np.einsum("ij, ij -> i", Z, Z)
    return np.sqrt(norms, norms)


def shared_borders(gdf, geo_field = "geometry"):
    rook = Rook.from_dataframe(gdf)
    weights = {
        idx: [
            gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length
            for nid in neighbours
        ]
        for idx, neighbours in rook.neighbors.items()
    }

    return W(rook.neighbors, weights).to_adjlist()


def knn_borders(gdf, geo_field="geometry"):
    rook = Rook.from_dataframe(gdf)
    weights = {
        idx: [
            gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length
            for nid in neighbours
        ]
        for idx, neighbours in rook.neighbors.items()
    }

    w_knn = KNN.from_dataframe(gdf, k=1)
    w = W(rook.neighbors, weights)
    w_attach = attach_islands(w, w_knn)

    return w_attach.to_adjlist()

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

    @classmethod
    def multi_to_single(cls, gdf, geo_field ="geometry"):
        gdf_singles = gdf[gdf[geo_field]].type == "Polygon"
        gdf_multis = gdf[gdf[geo_field]].type == "MultiPolygon"

        partial_seperation = partial(cls.__seperate, geo_field="geometry")
        sep = gdf_multis.apply(partial_seperation, axis=1).tolist()
        sep.append(gdf_singles)

        out = pd.concat(sep).reset_index(drop=True)
        out.crs = gdf.crs

        return out

    @staticmethod
    def __seperate(row, geo_field):
        df = pd.concat(
            [gpd.GeoDataFrame(row).T] * len(row[geo_field]), ignore_index=True
        )
        df[geo_field] = row[geo_field]
        return df

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


    def old_dorling(self, ratio=0.4, friction=0.25, iterations=99):
        def repel(x, row, xrepel, yrepel):
            if x["dist"] > 1.0:
                xrepel -= (
                    x["overlap"] * (x["geometry"].x - row["geometry"].x) / x["dist"]
                )
                yrepel -= (
                        x["overlap"] * (x["geometry"].y - row["geometry"].y) / x["dist"]
                )

            return xrepel, yrepel

        def attract(x, idx, borders, row, perimeter, xattract, yattract):
            # if (idx in borders["focal"].values) and (x.name in borders["neighbour"].values):
            #     x["overlap"] = (
            #         abs(x["overlap"]) * float(
            #             borders[(borders["focal"] == idx) & (borders["neighbor"] == x.name)]["weight"]
            #         ) / perimeter[idx]
            #     )

            # if (idx in borders["focal"].values) and (x.name in borders["neighbor"].values):
            if idx in borders["focal"].values:
                if x.name in borders[borders["focal"] == idx]["neighbor"].values:
                    x["overlap"] = (
                        abs(x["overlap"]) * float(
                            borders[(borders["focal"] == idx) & (borders["neighbor"] == x.name)]["weight"]
                        ) / perimeter[idx]
                    )

            xattract += x["overlap"] * (x["geometry"].x - row["geometry"].x) / x["dist"]
            yattract += x["overlap"] * (x["geometry"].y - row["geometry"].y) / x["dist"]

            return xattract, yattract

        borders = knn_borders(self.gdf)
        perimeter = self.gdf.length

        df = gpd.GeoDataFrame(
            self.gdf.drop(columns=self.geo_field), geometry=self.gdf.centroid
        )

        focal = np.stack(
            borders.merge(
                df[self.geo_field].map(np.array).to_frame(),
                left_on="focal",
                right_index=True,
            ).sort_index()[self.geo_field]
        )

        neighbor = np.stack(
            borders.merge(
                df[self.geo_field].map(np.array).to_frame(),
                left_on="neighbor",
                right_index=True,
            ).sort_index()[self.geo_field]
        )

        total_distance = np.sum(paired_distance(focal, neighbor))

        focal_radius = borders.merge(
            df[[self.value_field]],
            left_on="focal",
            right_index=True
        ).sort_index()[self.value_field]

        neighbor_radius = borders.merge(
            df[[self.value_field]],
            left_on="neighbor",
            right_index=True
        ).sort_index()[self.value_field]

        total_radius = np.sum(
            (focal_radius / np.pi) ** 0.5 + (neighbor_radius / np.pi) ** 0.5
        )

        scale = total_distance / total_radius

        df["radius"] = np.power(df[self.value_field] / np.pi, 0.5) * scale
        widest = df["radius"].max()

        for i in range(iterations):
            print(f"Iteration {i}")
            displacement = 0.0

            # For each geometry
            for idx, row in df.iterrows():
                xrepel = 0.0
                yrepel = 0.0
                xattract = 0.0
                yattract = 0.0
                closest = widest

                neighbours = df[
                    df.distance(row[self.geo_field]).between(
                        0, widest + row["radius"], inclusive='neither',
                    )
                ].copy()

                if len(neighbours) > 0:
                    neighbours["dist"] = neighbours[self.geo_field].distance(row[self.geo_field])

                    closest = widest if neighbours["dist"].min() > widest else neighbours["dist"].min()

                    neighbours["overlap"] = (neighbours["radius"] + row["radius"]) - neighbours["dist"]

                    for idy, rowy in neighbours.iterrows():
                        if rowy["overlap"] > 0.0:
                            xrepel, yrepel = repel(rowy, row, xrepel, yrepel)
                        else:
                            xattract, yattract = attract(rowy, idx, borders, row, perimeter, xattract, yattract)

                combined_attr = (xattract ** 2 + yattract ** 2) ** 0.5
                combined_repl = (xrepel ** 2 + yrepel ** 2) ** 0.5

                if combined_repl > closest:
                    xrepel = closest * xrepel / (combined_repl + 1.0)
                    yrepel = closest * yrepel / (combined_repl + 1.0)
                    combined_repl = closest

                if combined_repl > 0:
                    xtotal = (1.0 - ratio) * xrepel + ratio * (
                            combined_repl * xattract / (combined_attr + 1.0)
                    )
                    ytotal = (1.0 - ratio) * yrepel + ratio * (
                            combined_repl * yattract / (combined_attr + 1.0)
                    )
                else:
                    if combined_attr > closest:
                        xattract = closest * xattract / (combined_attr + 1.0)
                        yattract = closest * yattract / (combined_attr + 1.0)
                    xtotal = xattract
                    ytotal = yattract

                displacement += (xtotal ** 2 + ytotal ** 2) ** 0.5

                xvector = friction * xtotal
                yvector = friction * ytotal

                df.loc[idx, self.geo_field] = translate(
                    row[self.geo_field], xoff=xvector, yoff=yvector
                )

            displacement = displacement / len(df)

        return gpd.GeoDataFrame(
            data=df.drop(columns=["geometry", "radius"]),
            geometry=df.apply(lambda x: x["geometry"].buffer(x["radius"]), axis=1)
        )


    def dorling(self, ratio=0.4, friction=0.25, iterations=100):
        pass


    def diffusion(self):
        pass


    def fast_flow(self):
        pass
