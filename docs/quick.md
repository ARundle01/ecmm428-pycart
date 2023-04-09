# Quick Start

`pycart` is a Python library that allows the generation of Cartograms from a 
given GeoPandas GeoDataFrame. The library is designed around a single `Cartogram` 
class, with the generation algorithms being methods of said class.

## Tutorial
Creating a Cartogram begins with ensuring that your GeoDataFrame has the following 
basic structure:

```python
from pycart import cartogram
import geopandas as gpd

gdf = gpd.GeoDataFrame()

cart = cartogram.Cartogram(gdf, value_field="Population", id_field="Name", geometry_field="Geometry")
```