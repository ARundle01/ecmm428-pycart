# pycart - Python Cartogram Generation

`pycart` is a Python software package that allows generation of Cartograms from 
a GeoPandas GeoDataFrame.

The library provides a `cartogram` generator class, alongside individual methods 
based on popular Cartogram generation algorithms.

### Quick Install
`pycart` can be installed from PyPI, or the development version from TestPyPI.

Install from PyPI with:
```pycon
pip install pycart
```

Install the development version from TestPyPI:
```pycon
pip install -i https://test.pypi.org/simple/ pycart
```

For further details on installation, please see the [installation guide](install.md).

### Contributors

| Name                                        | Role           |
|---------------------------------------------|----------------|
| [Alex Rundle](https://github.com/ARundle01) | Lead Developer |

### License

This project is licensed under the BSD 3-Clause License. For further details, please see 
the [about section](about.md).

### Key Dependencies
To use `pycart`, you will need Python>=3.9, alongside the following 
fixed version dependencies:

| Package                                                    | Version      |
|------------------------------------------------------------|--------------|
| [Python](https://www.python.org/downloads/)                | >=3.9, <3.12 |
| [matplotlib](https://pypi.org/project/matplotlib/)         | 3.7.1        |
| [geopandas](https://pypi.org/project/geopandas/)           | 0.12.2       |
| [geojson](https://pypi.org/project/geojson/)               | 3.0.1        |
| [pandas](https://pypi.org/project/pandas/)                 | 1.5.3        |
| [numpy](https://pypi.org/project/numpy/)                   | 1.24.2       |
| [libpysal](https://pypi.org/project/libpysal/)             | 4.7.0        |
| [shapely](https://pypi.org/project/shapely/)               | 2.0.1        |
| [alive-progress](https://pypi.org/project/alive-progress/) | >3.1.1       |
