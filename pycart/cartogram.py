# Numerical Imports
import geopandas as gpd
import numpy as np

# Import shapely
from shapely import distance
from shapely.affinity import scale, translate

# Import progress bar for Dorling
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
    amount the neighbour overlaps with the focal region and $D_N$ is the distance
    from the neighbour to the focal region:

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
    # Get overlap, x and y differences and dist to neighbour
    overlap = neighbour['overlap']
    dx = neighbour['geometry'].x - focal['geometry'].x
    dy = neighbour['geometry'].y - focal['geometry'].y
    dist = neighbour['dist']

    # Subtract from repulsive forces
    xrepel -= overlap * dx / dist
    yrepel -= overlap * dy / dist

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
    # Create mask of if supplied focal and neighbour were originally neighbours
    mask = (borders["focal"] == idx) & (borders["neighbor"] == nb.name)

    # If focal and nb are neighbours
    if not borders[mask].empty:
        # Scale overlap proportional to border weight
        nb["overlap"] = (np.abs(nb["overlap"]) * float(borders[mask]["weight"]) / perimeter[idx])

    # Get overlap, x and y differences and dist to neighbour
    overlap = nb["overlap"]
    dx = nb["geometry"].x - focal["geometry"].x
    dy = nb["geometry"].y - focal["geometry"].y
    dist = nb["dist"]

    # Add to attractive forces
    xattract += overlap * dx / dist
    yattract += overlap * dy / dist

    return xattract, yattract


class Cartogram:
    def __init__(self, gdf, value_field, id_field, geometry_field='geometry'):
        """
        A Cartogram object acts as a generator on a specific dataset, through which
        generation algorithms can be run.

        ###**Parameters**

        - **gdf  :  *[geopandas.GeoDataFrame](https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.html)*** - The dataset that you want to apply cartogram techniques to
        - **value_field  :  *String*** - Field in the dataset to apply cartogram techniques to
        - **id_field  :  *String*** - Field of the identifier for each region in the dataset
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
        self.id_field = id_field

    def non_contiguous(self, size_value=1.0):
        """
        Calculates and returns a Non-Contiguous cartogram in the form of a GeoDataFrame.

        A Non-Contiguous cartogram [[1]](cartogram.md#References) is created by scaling each
        region in-place by a specific density (value field divided by geographical area), about an
        anchor region. The anchor region is usually the region with the highest density.

        Suppose we have an anchor unit $H$ of area $A_H$ and value $V_H$ and a miscellaneous
        region $J$, with area $A_J$ and value $V_J$. The scaling value applied to $J$ is:

        $(\sqrt{V_H / A_H})^{-1} \cdot \sqrt{V_J / A_J}$

        ###**Parameters**

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
        non_con = cart.non_contiguous(size_value=1.0)

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
        # Create copy of supplied gdf to prevent permanent data changes
        geodf = self.gdf[[self.value_field, self.id_field, self.geo_field]].copy()

        # Get centroids of all regions
        geodf["centroid"] = geodf[self.geo_field].centroid

        # Calculate value-density
        geodf["density"] = geodf[self.value_field] / geodf.area

        # Anchor is the region with the highest density
        anchor = geodf.loc[geodf["density"].idxmax(), "density"]

        # Scale = \sqrt(density / anchor density)
        # size_value is a simple multiplier to enhance scaling if outputs are weak
        geodf["scale"] = (geodf["density"] / anchor).pow(0.5) * size_value

        # Scale geometries by scale value about centroid
        new_geo = []
        for _, geo in geodf.iterrows():
            scaled = scale(geo[self.geo_field], geo['scale'], geo['scale'], origin=geo['centroid'])
            new_geo.append(scaled)

        geodf = geodf.drop(columns=['density', 'centroid'])

        return gpd.GeoDataFrame(geodf, geometry=new_geo)

    def dorling(self, iterations=100, ratio=0.4, friction=0.5, stop=None):
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
        # Create copy of supplied gdf to prevent permanent data changes
        geodf = self.gdf[[self.value_field, self.id_field, self.geo_field]].copy()

        # Calculate borders
        borders, _ = border.get_borders(geodf)

        perimeter = geodf.length

        # Replace geometries with centroids
        regions = gpd.GeoDataFrame(geodf.drop(columns=self.geo_field), geometry=geodf.centroid)

        # Get focal regions and neighbour regions
        focals = regions[self.geo_field].loc[borders["focal"]].values
        neighbours = regions[self.geo_field].loc[borders["neighbor"]].values

        # Calculate sum of paired_distances between all regions and neighbours
        total_distance = np.sum(_paired_distance(focals, neighbours))

        # Calculate area of each region, where area is equal to the value_field
        focal_area = borders.merge(regions[[self.value_field]], left_on="focal", right_index=True).sort_index()[self.value_field]
        neighbour_area = borders.merge(regions[[self.value_field]], left_on="neighbor", right_index=True).sort_index()[self.value_field]

        # Calculate the sum of unscaled radii for each region
        # The radius for a region is found using r = sqrt(area / pi)
        total_radius = np.sum(np.sqrt(focal_area / np.pi) + np.sqrt(neighbour_area / np.pi))

        # Calculate scale coefficient
        scale = total_distance / total_radius

        # Calculate scaled radius by multiplying original radius by scale
        regions["radius"] = np.power(regions[self.value_field] / np.pi, 0.5) * scale

        # Get widest region
        widest = regions["radius"].max()

        with alive_bar(iterations, title='Making Dorling Cartogram') as bar:
            for i in range(iterations):
                # print(f"Starting Iteration: {i}")
                bar()

                if stop is not None:
                    if i == stop:
                        break

                for idx, region in regions.iterrows():
                    xrepel = yrepel = xattract = yattract = 0.0
                    closest = widest

                    # Get all regions that are within a radius r, where r is the widest radius + current region radius
                    neighbours = regions[regions.distance(region[self.geo_field]).between(0, widest + region["radius"], inclusive="neither")].copy()

                    if len(neighbours) > 0:
                        # Calculate distance between each neighbour and the current region
                        neighbours["dist"] = neighbours[self.geo_field].distance(region[self.geo_field])

                        # Get the closest distance between the current region and all neighbours
                        closest = widest if neighbours["dist"].min() > widest else neighbours["dist"].min()

                        # Calculate overlap for all neighbours
                        neighbours["overlap"] = (neighbours["radius"] + region["radius"]) - neighbours["dist"]

                        # Calculate repulsive and attractive forces
                        for idy, nb in neighbours.iterrows():
                            if nb["overlap"] > 0.0:
                                xrepel, yrepel = _repel(nb, region, xrepel, yrepel)
                            else:
                                xattract, yattract = _attract(nb, borders, idx, region, perimeter, xattract, yattract)

                    # Calculate the distance covered by attractive and repulsive forces
                    attract_dist = np.hypot(xattract, yattract)
                    repel_dist = np.hypot(xrepel, yrepel)

                    # If region is farther away than the closest region, scale vector
                    if repel_dist > closest:
                        xrepel, yrepel = closest * np.array((xrepel, yrepel)) / (repel_dist + 1.0)
                        repel_dist = closest

                    # Calculate new position of region based on forces
                    if repel_dist > 0:
                        xtotal = (1 - ratio) * xrepel + ratio * (
                                repel_dist * xattract / (attract_dist + 1.0))
                        ytotal = (1 - ratio) * yrepel + ratio * (
                                repel_dist * yattract / (attract_dist + 1.0))
                    else:
                        if attract_dist > closest:
                            xattract, yattract = closest * np.array((xattract, yattract)) / (attract_dist + 1.0)
                        xtotal, ytotal = xattract, yattract

                    # Calculate velocity vector by applying friction to total force vectors
                    xvector, yvector = friction * np.array((xtotal, ytotal))

                    # Update current region position based on calculate xy vectors
                    regions.loc[idx, self.geo_field] = translate(
                        region[self.geo_field], xoff=xvector, yoff=yvector
                    )

                    # Create the circle regions by adding a buffer around the centroid of
                    # region using the calculated radius
                    buffered_geos = []
                    for _, row in regions.iterrows():
                        buffered_geo = row["geometry"].buffer(row['radius'])
                        buffered_geos.append(buffered_geo)

        return gpd.GeoDataFrame(
            data=regions.drop(columns=["geometry", "radius"]),
            geometry=buffered_geos
        )
