from libpysal.weights import Queen, W


def get_borders(gdf, geo_field="geometry"):
    """
    Calculates bordering regions for all regions using Queen contiguity (i.e. neighbours share
    either a corner or an edge).

    :param gdf: GeoDataFrame containing regions
    :param geo_field: Field within gdf that contains the geometric data
    :return: Border Adjacency list and list of Islands (disconnected components)
    """
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

        borders = W(wq.neighbors, weights, silence_warnings=True).to_adjlist()

        borders = borders.astype({'focal': int, 'neighbor': int})

        return borders, islands
    else:
        weights = {
            idx: [
                gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length
                for nid in neighbors
            ]
            for idx, neighbors in wq.neighbors.items()
        }

        borders = W(wq.neighbors, weights).to_adjlist()

        borders = borders.astype({'focal': int, 'neighbor': int})

        return borders, islands
