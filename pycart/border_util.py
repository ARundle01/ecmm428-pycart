from libpysal.weights import Queen, W
import pandas as pd


def get_borders(gdf, geo_field="geometry"):
    """
    Calculates bordering regions for all regions using Queen contiguity (i.e. neighbours share
    either a corner or an edge).

    Similar to a [Moore neighbourhood](https://en.wikipedia.org/wiki/Moore_neighborhood), Queen
    contiguity determines a neighbour as a region that shares an edge or vertex with the focal region.

    Weights are assigned to each focal, neighbour pair using [`libpysal.weights.Queen`](https://pysal.org/libpysal/generated/libpysal.weights.Queen.html).
    The weights are based on the amount of a focal region's perimeter is occupied by the neighbour region.

    ###**Parameters**

    - **gdf  :  *[geopandas.GeoDataFrame](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html)*** - A dataset that you want to find all borders for
    - **geo_field  :  *string, optional, default 'geometry'*** - Field containing the geometries of each region

    ###**Returns**

    - ***[pandas.DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)*** - DataFrame where each row contains
    the focal node, its neighbour and the Queen contiguity weight.

    """
    islands = None
    wq = Queen.from_dataframe(gdf, silence_warnings=True)

    # If there are disconnected regions (islands)
    if len(wq.islands) > 0:
        # Remove the islands
        islands = wq.islands
        temp_gdf = gdf.drop(islands, axis=0)

        # Reset the index after removing the islands
        WQ = wq.to_adjlist(drop_islands=True)
        WQ.reset_index(inplace=True, drop=True)
        wq = W.from_adjlist(WQ)

        # Compute the length of intersection between current region and neighbours
        weights = {
            idx: [
                temp_gdf.loc[idx, geo_field].intersection(temp_gdf.loc[nid, geo_field]).length
                for nid in neighbors
            ]
            for idx, neighbors in wq.neighbors.items()
        }

        # Create the borders DataFrame
        borders = pd.DataFrame(columns=['focal', 'neighbor', 'weight'])

        # For each focal region, extract weights for all neighbours of focal
        for focal, neighbour_list in wq.neighbors.items():
            weight_list = weights[focal]

            temp_df = pd.DataFrame({'focal': focal, 'neighbor': neighbour_list, 'weight': weight_list})

            # Add to DataFrame
            borders = pd.concat([borders, temp_df], ignore_index=True)

        # Convert region numbers to ints
        borders = borders.astype({'focal': int, 'neighbor': int})

        return borders, islands
    else:
        # Compute the length of intersection between current region and neighbours
        weights = {
            idx: [
                gdf.loc[idx, geo_field].intersection(gdf.loc[nid, geo_field]).length
                for nid in neighbors
            ]
            for idx, neighbors in wq.neighbors.items()
        }

        # Create the borders DataFrame
        borders = pd.DataFrame(columns=['focal', 'neighbor', 'weight'])

        # For each focal region, extract weights for all neighbours of focal
        for focal, neighbour_list in wq.neighbors.items():
            weight_list = weights[focal]

            temp_df = pd.DataFrame({'focal': focal, 'neighbor': neighbour_list, 'weight': weight_list})

            # Add to DataFrame
            borders = pd.concat([borders, temp_df], ignore_index=True)

        # Convert region numbers to ints
        borders = borders.astype({'focal': int, 'neighbor': int})

        return borders, islands
