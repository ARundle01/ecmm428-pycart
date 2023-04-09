# Third Party Imports
import geopandas as gpd
import numpy as np

from shapely import distance
from shapely.affinity import scale, translate

# Project Library
import pycart.border_util as border


def _paired_distance(X, Y):
    """
    Calculates the pairwise Euclidean distance between two arrays of Points
    :param X: ndarray of Points
    :param Y: ndarray of Points
    :return: ndarray of Euclidean distances
    """
    distances = np.array([distance(x, y) for x, y in zip(X, Y)])
    return distances


class Cartogram:
    def __init__(self, gdf, value_field, id_field=None, geometry_field='geometry'):
        """
        Initialise the Cartogram generator for a given dataset.

        :param gdf: Dataset that you want to transform
        :param value_field: Field in the dataframe to apply cartogram techniques to, e.g. Population
        :param id_field: Field of any identifier for each region in gdf
        :param geometry_field: Field containing geometries of each region
        """
        # Initialise general GeoDataFrame attributes
        self.gdf = gdf
        self.value_field = value_field
        self.geo_field = geometry_field

        if not id_field:
            self.gdf['id_field'] = self.gdf.index
            self.id_field = "id_field"
        else:
            self.id_field = id_field

    def non_contiguous(self, position="centroid", size_value=1.0):
        """
        Calculates and returns a Non-Contiguous cartogram in the form of a GeoDataFrame.

        :param position: The position to use when scaling each region
        :param size_value: Simple multiplier to the scaling of each region
        :return: GeoDataFrame representation of the non-contiguous cartogram
        """
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

        geodf["scale"] = (1.0 / np.power(anchor, 0.5)) * np.power(geodf[self.value_field] / geodf.area,
                                                                  0.5) * size_value

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
        """
        Runs the Dorling cartogram algorithm and returns the generated cartogram in the form of a GeoDataFrame.

        :param ratio: Ratio of Attraction force to Repulsion force; the larger the ratio, the more attraction force
        :param friction: Determines the strength of x and y vectors; acts as a simple multiplier to x and y vectors
        :param iterations: The number of iterations to run the Dorling algorithm for.
        :param stop: A given iteration to halt computation at.
        :return: GeoDataFrame representation of a Dorling cartogram
        """
        def repel(neighbour, focal, xrepel, yrepel):
            """
            Calculates the repulsive forces in X and Y components.

            :param neighbour: List of Neighbour points
            :param focal: List of Focal points
            :param xrepel: existing repulsive force in X direction
            :param yrepel: existing repulsive force in Y direction
            :return: New X and Y direction repulsive forces
            """
            xrepel -= (
                    neighbour["overlap"] * (neighbour["geometry"].x - focal["geometry"].x) / neighbour["dist"]
            )
            yrepel -= (
                    neighbour["overlap"] * (neighbour["geometry"].y - focal["geometry"].y) / neighbour["dist"]
            )

            return xrepel, yrepel

        def attract(nb, borders, idx, focal, perimeter, xattract, yattract):
            """
            Calculates and updates the attractive forces in X and Y components

            :param nb: An individual Neighbour
            :param borders: List of borders for all regions in self.gdf
            :param idx: The index of the current region
            :param focal: The current region
            :param perimeter: The perimeter of the current region
            :param xattract: Existing attractive forces in the X direction
            :param yattract: Existing attractive forces in the Y direction
            :return: New X and Y direction attractive forces
            """
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

        borders, islands = border.get_borders(self.gdf)
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

        total_distance = np.sum(_paired_distance(focal, neighbour))

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

            if stop is not None:
                if i == stop:
                    break

            for idx, region in regions.iterrows():
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

        return gpd.GeoDataFrame(
            data=regions.drop(columns=["geometry", "radius"]),
            geometry=regions.apply(lambda x: x["geometry"].buffer(x["radius"]), axis=1)
        )
