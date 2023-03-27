from libpysal.weights import Queen, W


def get_borders(gdf, geo_field="geometry"):
    islands = None
    wq = Queen.from_dataframe(gdf, silence_warnings=True)

    if len(wq.islands) > 0:
        islands = wq.islands

        temp_gdf = gdf.drop(islands, axis=0)

        WQ = wq.to_adjlist(drop_islands=True)
        WQ.reset_index(inplace=True, drop=True)
        wq = W.from_adjlist(WQ)

        weights = {
            idx: [
                temp_gdf.loc[idx, geo_field].intersection(temp_gdf.loc[nid, geo_field]).length
                for nid in neighbors
            ]
            for idx, neighbors in wq.neighbors.items()
        }

        return W(wq.neighbors, weights, silence_warnings=True).to_adjlist(), islands
    else:
        weights = {
            idx: [
                gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length
                for nid in neighbors
            ]
            for idx, neighbors in wq.neighbors.items()
        }

        return W(wq.neighbors, weights).to_adjlist(), islands
