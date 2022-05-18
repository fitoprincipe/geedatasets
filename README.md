### Google Earth Engine Datasets
(_this package was a module in [geetools](https://github.com/gee-community/gee_tools) 
until vesion 0.5. Now it is a package by itself but still depends on geetools_)

In this package you can find hardcoded information about Google Earth Engine 
available datasets, such as Landsat and Sentinel, and also methods for 
manipulating those datasets, such as cloud masking and indexes calculations.

### WARNING
New version 0.2.0 incorporates Landsat Collection 2 and makes it the default
Landsat collection. So this could be a braking change

## Install
> pip install geedatasets

## Usage
All available datasets are registerd in the global variable `ALL` but there are
two getter methods:

### `fromId`
This method gets the dataset from its `id` if present in `ALL`

```python
import geedatasets
# Sentinel 2 TOA
s2toa = geedatasets.fromId('COPERNICUS/S2')
```

### `fromShortName`
This method gets the dataset from its `short_name` if present in `ALL`

```python
import geedatasets
# Landsat 8 Surface Reflectance (Collection 2 since version 0.2.0)
l8sr = geedatasets.fromShortName('L8SR')
```

## Available Datasets
- Landsat
- Sentinel 2
- Sentinel 1