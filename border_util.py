import geopandas as gpd
import numpy as np
import pandas as pd

from libpysal.weights import Queen, W, Kernel


def bad_borders(gdf, geo_field="geometry"):
    wq = Queen.from_dataframe(gdf)

    if len(wq.islands) > 0:
        wk = Kernel.from_dataframe(gdf, fixed=False, k=5)
        WK  = wk.to_adjlist()
        WK = WK.loc[WK["focal"].isin(wq.islands)]
        WK = WK[WK["focal"] != WK["neighbor"]]
        WK = WK[WK["weight"] > 0.001]

        WQ = wq.to_adjlist(drop_islands=True)

        COMBINED = pd.concat([WQ, WK]).sort_values(by=['focal', 'neighbor']).reset_index(drop=True)

        combined_W = W.from_adjlist(COMBINED)

        weights = {
            idx: [
                gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length
                for nid in neighbors
            ]
            for idx, neighbors in combined_W.neighbors.items()
        }

        for idx, neighbors in combined_W.neighbors.items():
            idx_weights = []
            for nid in neighbors:
                # if gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length == 0:
                #     print(f"{idx}: {nid}")
                print(gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]))

        # print(len(combined_W.neighbors.items()))

        # TODO: to_adjlist() throws ValueError: All arrays must be of the same length
        return W(combined_W.neighbors, weights).to_adjlist()
    else:
        weights = {
            idx: [
                gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length
                for nid in neighbors
            ]
            for idx, neighbors in wq.neighbors.items()
        }

        return W(wq.neighbors, weights).to_adjlist()

    # weights = {
    #     idx: [
    #         gdf.loc[idx, geo_field].intersection(gdf.loc(nid, geo_field)).length
    #         for nid in neighbors
    #     ]
    #     for idx, neighbors in queen.neighbors.items()
    # }
    #
    # return W(queen.neighbors, weights).to_adjlist()


def get_borders(gdf, geo_field="geometry"):
    islands = None
    wq = Queen.from_dataframe(gdf)

    # if len(wq.islands) > 0:
    #     islands = wq.islands
    #
    #     temp_gdf = gdf.drop(islands, axis=0)
    #
    #     WQ = wq.to_adjlist(drop_islands=True)
    #     WQ.reset_index(inplace=True, drop=True)
    #     wq = W.from_adjlist(WQ)
    #
    #     weights = {
    #         idx: [
    #             temp_gdf.loc[idx, geo_field].intersection(temp_gdf.loc[nid, geo_field]).length
    #             for nid in neighbors
    #         ]
    #         for idx, neighbors in wq.neighbors.items()
    #     }
    # else:
    weights = {
        idx: [
            gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length
            for nid in neighbors
        ]
        for idx, neighbors in wq.neighbors.items()
    }

    return W(wq.neighbors, weights).to_adjlist(), islands
