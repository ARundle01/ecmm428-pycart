# Third Party Imports
import geopandas as gpd
import numpy as np

from shapely import distance
from shapely.affinity import scale, translate

from alive_progress import alive_bar

# Project Library
import pycart.border_util as border


def _paired_distance(X, Y):
    """
    Calculates the pairwise Euclidean distance between two *array_like* of [`shapely.Point`](https://shapely.readthedocs.io/en/stable/reference/shapely.Point.html),
    using [`shapely.distance()`](https://shapely.readthedocs.io/en/stable/reference/shapely.distance.html).

    Both **X** and **Y** should be the same length.

    ###**Parameters**

    - **X, Y  :  *array_like of `shapely.Point`*** - The two arrays to calculate distance between.

    ###**Returns**

    - ***numpy.ndarray*** - array of `float` distances
    """
    distances = np.array([distance(x, y) for x, y in zip(X, Y)])
    return distances


def _repel(neighbour, focal, xrepel, yrepel):
    """
    Calculates and updates the repulsive force being applied onto a given region, from a
    given neighbour of said region.

    The repulsive force $F$ is split into $x$ and $y$ components, $F_x$ and $F_y$
    respectively. These forces are calculated as the following, where $O$ is the
    amount the neighbour overlaps with the focal region:

    $F_x = O \cdot (N_x - M_x) / D_N$

    $F_y = O \cdot (N_y - M_y) / D_N$

    The supplied $x$ and $y$ forces are updated by subtracting the new forces.

    ###**Parameters**

    - **neighbour  :  *[pandas.Series](https://pandas.pydata.org/docs/reference/api/pandas.Series.html)*** - A given Neighbour of **focal**
    - **focal  :  *[pandas.Series](https://pandas.pydata.org/docs/reference/api/pandas.Series.html)*** - The current, focal region
    - **xrepel  :  *float*** - The current repulsive force in the $x$ direction
    - **yrepel  :  *float*** - The current repulsive force in the $y$ direction

    ###**Returns**

    - **xrepel  :  *float*** - The new repulsive force in the $x$ direction.
    - **yrepel  : *float*** - The new repulsive force in the $y$ direction.
    """
    xrepel -= (
            neighbour["overlap"] * (neighbour["geometry"].x - focal["geometry"].x) / neighbour["dist"]
    )
    yrepel -= (
            neighbour["overlap"] * (neighbour["geometry"].y - focal["geometry"].y) / neighbour["dist"]
    )

    return xrepel, yrepel


def _attract(nb, borders, idx, focal, perimeter, xattract, yattract):
    """
    Calculates and updates the attractive force being applied to a given region, towards a given
    neighbour region.

    Before the attractive forces are calculated, the overlap $O$ amount for a neighbour is scaled as
    such, where $W_{FN}$ is the Queen contiguity weight and $P_F$ is the perimeter of the focal region:

    $O_{new} = (| O_{original} | * W_{NF}) / P_F$

    The attractive force $A$ is split into $x$ and $y$ components, $A_x$ and $A_y$
    respectively. These forces are calculated as the following, where $O$ is the
    amount the neighbour overlaps with the focal region:

    $A_x = O_x \cdot (N_x - M_x) / D_N$

    $A_y = O_y \cdot (N_y - M_y) / D_N$

    The supplied $x$ and $y$ forces are updated by adding the new forces.

    ###**Parameters**

    - **nb  : *[pandas.Series](https://pandas.pydata.org/docs/reference/api/pandas.Series.html)*** - A given Neighbour of **focal**
    - **borders  :  *[pandas.DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)*** - DataFrame of all border weights
    - **idx  :  *int*** - The index of the current region
    - **focal  :  *[pandas.Series](https://pandas.pydata.org/docs/reference/api/pandas.Series.html)*** - The current, focal region
    - **perimeter :  *[pandas.Series](https://pandas.pydata.org/docs/reference/api/pandas.Series.html)*** - The perimeters of all regions
    - **xattract  :  *float*** - The current attractive force in the $x$ direction
    - **yattract  :  *float*** - The current attractive force in the $y$ direction

    ###**Returns**

    - **xattract  :  *float*** - The new attractive force in the $x$ direction.
    - **yattract  : *float*** - The new attractive force in the $y$ direction.
    """
    if sum((borders["focal"] == idx) & (borders["neighbor"] == nb.name)) == 1:
        nb["overlap"] = (
                np.abs(nb["overlap"])
                * float(borders[(borders["focal"] == idx) & (borders["neighbor"] == nb.name)]["weight"])
                / perimeter[idx]
        )

    xattract += nb["overlap"] * (nb["geometry"].x - focal["geometry"].x) / nb["dist"]
    yattract += nb["overlap"] * (nb["geometry"].y - focal["geometry"].y) / nb["dist"]

    return xattract, yattract


class Cartogram:
    def __init__(self, gdf, value_field, id_field=None, geometry_field='geometry'):
        """
        A Cartogram object acts as a generator on a specific dataset, through which
        generation algorithms can be run.

        ###**Parameters**

        - **gdf  :  *[geopandas.GeoDataFrame](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html)*** - The dataset that you want to apply cartogram techniques to
        - **value_field  :  *String*** - Field in the dataset to apply cartogram techniques to
        - **id_field  :  *String, optional, default None*** - Field of the identifier for each region in the dataset
        - **geometry_field  :  *String, optional, default 'geometry'*** - Field containing the geometries of each region

        ###**Example**
        ```python
        from pycart import cartogram
        import geopandas as gpd

        my_geodf = gpd.read_file("path/to/dataset.csv")

        cart = cartogram.Cartogram(my_geodf, value_field='Population', id_field='Name', geo_field='Geometry')
        ```

        ---
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

        A Non-Contiguous cartogram [[1]](cartogram.md#References) is created by scaling each
        region in-place by a specific density (value field divided by geographical area), about an
        anchor region. The anchor region is usually the region with the highest density.

        Suppose we have an anchor unit $H$ of area $A_H$ and value $V_H$ and a miscellaneous
        region $J$, with area $A_J$ and value $V_J$. The scaling value applied to $J$ is:

        $\sqrt{(A_H / A_J) \cdot (V_J / V_H)}$

        ###**Parameters**

        - **position  :  *{'centroid', 'centre'}, optional, default 'centroid'*** - Apply scaling based on the region's centroid or the centre of the entire map.
        - **size_value  :  *int, optional, default 1.0*** - A simple multiplier to scaling; larger values accentuate scaling.

        ###**Returns**

        - ***[geopandas.GeoDataFrame](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html)*** - Scaled form of original dataset

        ###**Example**
        ```python
        from pycart import cartogram
        import matplotlib.pyplot as plt
        import geopandas as gpd

        # Load dataset into Cartogram generator
        my_geodf = gpd.read_file("path/to/dataset.csv")
        cart = cartogram.Cartogram(my_geodf, value_field='Population', id_field='Name', geo_field='Geometry')

        # Create Non-Contiguous cartogram
        non_con = cart.non_contiguous(position='centroid', size_value=1.0)

        # Plot data
        fig, ax = plt.subplots(1, figsize=(4, 4))
        ax.axis('equal')
        ax.axis('off')

        my_geodf.plot(color='w', ax=ax, alpha=0.8, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)
        non_con.plot(color='r', ax=ax, edgecolor='0', linewidth=0.1, legend=False)

        plt.show()
        ```

        ---
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

    def dorling(self, iterations=100, ratio=0.4, friction=0.25, stop=None):
        """
        Runs the Dorling cartogram algorithm and returns the generated cartogram in the form of a GeoDataFrame.

        A Dorling cartogram [[2]](cartogram.md#References) represents each region as a circle, with the
        radius being proportional to the density of the given region. The regions are then moved via a
        gravity-like force model.

        The radius of a given region $J$ is calculated over multiple steps, with Steps 1 and 2 being carried out by
        [`border_util.get_borders()`](border.md):

        1. Use Queen contiguity to find all neighbours of a region $J$; Queen contiguity acts similarly to
        a [Moore neighbourhood](https://en.wikipedia.org/wiki/Moore_neighborhood), where a neighbour is a region
        that shares an edge or a vertex.
        2. Assign weights to each neighbour of $J$ based on the length of perimeter that the neighbour occupies.
        3. Calculate the sum $D$, of the pairwise distance (see [`_paired_distance`](cartogram.md#Helpers)) from a region to all of its neighbours, for all regions.
        4. Calculate the sum $R$, of radii as:
        $R = \sum (\sqrt{V_M / \pi} + \sqrt{V_N / \pi})$ for all $N$ neighbours and for all $M$ regions

        5. Calculate the radius scaling coefficient $k$, as:
        $k = D / R$

        6. Calculate the radius of region $J$:
        $r_J = k \sqrt{V_J / \pi}$

        See the helper functions [`_attract()`](cartogram.md#Helpers) and [`_repel()`](cartogram.md#Helpers) for more details on how the forces are
        calculated.

        ###**Parameters**

        - **iterations  :  *int, default 100*** - The number of iterations to run the Dorling algorithm for.
        - **ratio  : *float, optional, default 0.4*** - Ratio of attractive force to repulsive force; the larger the ratio, the more attractive force is applied.
        - **friction  :  *float, optional, default 0.25*** - Determines the strength of $x$ and $y$ vectors; acts as a simple multiplier to $x$ and $y$ vectors.
        - **stop  :  *int, optional, default None*** - A given iteration at which to halt computation and return the cartogram

        ###**Returns**

        - ***[geopandas.GeoDataFrame](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html)*** - Scaled and Distributed circles based on original dataset

        ###**Example**
        ```python
        from pycart import cartogram
        import matplotlib.pyplot as plt
        import geopandas as gpd

        # Load dataset into Cartogram generator
        my_geodf = gpd.read_file("path/to/dataset.csv")
        cart = cartogram.Cartogram(my_geodf, value_field='Population', id_field='Name', geo_field='Geometry')

        # Run Dorling algorithm for 100 iterations
        dorling = cart.dorling(iterations=100, stop=None)

        # Plot data
        fig, ax = plt.subplots(1, figsize=(4, 4))
        ax.axis('equal')
        ax.axis('off')

        my_geodf.plot(color='w', ax=ax, alpha=0.8, zorder=0,  edgecolor='0', linewidth=0.1, legend=False)
        dorling.plot(color='w', ax=ax, alpha=0.8, zorder=0, edgecolor='0', linewidth=0.1, legend=False)

        plt.show()
        ```
        """

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

        with alive_bar(iterations, title='Making Dorling Cartogram') as bar:
            for i in range(iterations):
                # print(f"Starting Iteration: {i}")
                bar()
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
                                xrepel, yrepel = _repel(nb, region, xrepel, yrepel)
                            else:
                                xattract, yattract = _attract(nb, borders, idx, region, perimeter, xattract, yattract)

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
