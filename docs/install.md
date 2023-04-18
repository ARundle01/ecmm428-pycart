# Installing `pycart`

Before installing `pycart`, make sure you have the latest version of pip:
```pycon
py -m pip install --upgrade pip
```

## Installing from PyPI
The latest stable version of `pycart` can be installed from [PyPI](https://pypi.org/project/pycart/):

```pycon
pip install pycart
```

## Installing from TestPyPI
The most recent development version of `pycart` can be installed from [TestPyPI](https://test.pypi.org/project/pycart/):

```pycon
pip install -i https://test.pypi.org/simple/ pycart
```

## Installing from Source
Download the source (tar.gz) from [PyPI](https://pypi.org/project/pycart/) or from the latest development version on 
[GitHub](https://github.com/ARundle01/ecmm428-pycart/tree/dev).

Unpack and change to the source directory (which should contain the README.md); for the purposes of this guide, 
the source folder will be referred to by the filepath `./pycart`.

`pycart` was originally built and packaged using [Poetry](https://python-poetry.org/), so `pycart` can be installed 
from source using either `poetry` or regularly using `python`.

### Installing with `poetry`
First, make sure you have the latest version of Poetry installed (please see their 
[installation guide](https://python-poetry.org/docs/)) and ensure that you have configured Poetry 
to create Virtual Environments in the working directory:

```shell
poetry config virtualenvs.in-project true
```

Navigate to the extracted source folder:

```
cd './pycart'
```

Create a new virtual environment using Poetry:
```shell
poetry env use python
```
This will create a new environment in the pycart directory that uses the latest version of Python.

Now, you can install all dependencies from source and subsequently install `pycart`:

```shell
poetry install
```

### Installing with `pip`
If you do not have Poetry installed, you can install `pycart` from source using `pip`.

First, make sure you have the latest version of `pip` installed:
```pycon
py -m pip install --upgrade pip
```

Next, create a new virtual environment in the source directory:
```pycon
python -m venv .venv
```

Activate the virtual environment using the following command:
```pycon
.venv/Scripts/activate.bat
```

Install the dependencies required by `pycart` first:
```pycon
py -m pip install -r requirements.txt
```

Finally, install the `pycart` package from the source directory:
```pycon
pip install -e './pycart'
```

## Key Dependencies
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

