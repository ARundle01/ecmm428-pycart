ecmm428 - The `pycart` package
==============================

![GitHub](https://img.shields.io/github/license/ARundle01/ecmm428-pycart)
![GitHub code size in bytes](https://img.shields.io/github/repo-size/ARundle01/ecmm428-pycart)
![GitHub last commit](https://img.shields.io/github/last-commit/ARundle01/ecmm428-pycart)
![GitHub top language](https://img.shields.io/github/languages/top/ARundle01/ecmm428-pycart)
[![Documentation Status](https://readthedocs.org/projects/ecmm428-pycart/badge/?version=latest)](https://ecmm428-pycart.readthedocs.io/en/latest/?badge=latest)
![PyPI](https://img.shields.io/pypi/v/pycart)

## Contributors
Contributors to this project are
- Alex Rundle (ARundle01)

## License
This project is licensed under the BSD 3-Clause License

SPDX-License-Identifier: BSD-3-Clause

## Installation
### Using Pip
Make sure you have the latest version of PIP installed:
```pycon
py -m pip install --upgrade pip
```

This package is currently hosted on PyPi and TestPyPi, as it is currently in unstable development. 

To install the latest stable version
and use this package from PyPi:
```pycon
py -m pip install pycart
```

To install
and use the latest development version from TestPyPi:
```pycon
py -m pip install -i https://test.pypi.org/simple/ pycart
```

## Pre-requisites
All code in this Repository and Package was developed and tested on Windows 10.

This Repository requires:
- **Python >=3.9.7, <3.10**

## Read the Docs!
Further information about `pycart` can be found [here](https://ecmm428-pycart.readthedocs.io/en/latest/), 
including further installation and dependency details, a quickstart guide and API reference.

You can also read the research paper that accompanies this package [here](./ecmm428-pycart-a-library-for-cartogram-generation-in-python-with-acknowledgements.pdf).

### Key Dependencies
To use `pycart`, you will need Python>=3.9, alongside the following 
fixed version dependencies:

| Package                                                    | Version      |
|------------------------------------------------------------|--------------|
| [Python](https://www.python.org/downloads/)                | >=3.9, <3.10 |
| [matplotlib](https://pypi.org/project/matplotlib/)         | 3.7.1        |
| [geopandas](https://pypi.org/project/geopandas/)           | 0.12.2       |
| [geojson](https://pypi.org/project/geojson/)               | 3.0.1        |
| [pandas](https://pypi.org/project/pandas/)                 | 1.5.3        |
| [numpy](https://pypi.org/project/numpy/)                   | 1.24.2       |
| [libpysal](https://pypi.org/project/libpysal/)             | 4.7.0        |
| [shapely](https://pypi.org/project/shapely/)               | 2.0.1        |
| [alive-progress](https://pypi.org/project/alive-progress/) | >3.1.1       |