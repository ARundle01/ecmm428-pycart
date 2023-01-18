import geojson
import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from shapely.affinity import scale, translate
from shapely.geometry import Polygon
from libpysal.weights.contiguity import Rook, W


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

    def non_contiguous(self, position="centroid", anchor_rank=1, anchor_id=None):
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

        # for g in self.gdf.iterrows():
        #     print(g[1]["Name"])

        geodf["scale"] = (1.0 / np.power(anchor, 0.5)) * np.power(geodf[self.value_field] / geodf.area, 0.5)

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
