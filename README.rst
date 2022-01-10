==================
jet-stream-metrics
==================
To-Do:

- ADD cf_xarray (see: https://cf-xarray.readthedocs.io/en/latest/index.html)
- THINK ABOUT way of seperating main.py/experiments with jsmetrics
- ADD check to see if data has been loaded into memory (versus dask)
- way of dealing with data from different sources (some sort of data translator module or maybe included in tests)
  - for example what if 'v' or 'v-wind' is passed to func instead of 'va'
  - for example what if 'mbar' or 'model levels' instead of 'plev'
  - zg and z do not line up because of different units (geopotential_by_the_Earths_gravitational_acceleration = 9.80665, https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height)
- ADD way of subsetting longitude if it wraps around 0-360
- IMPROVE readme(s) and add links (DOI?)m (https://github.com/othneildrew/Best-README-Template/blob/master/README.md)
- ADD lat and lon subsetting to jetstream_metrics_dict
- IMPROVE performance. Look into bash arrays and running multiple metrics
- ADD benchmark speed tests for each metric based on standard sizes of the dataset
- argparse to parse a list for additional arguments (do more user testing)
- Check that methods using just U account for negative values (maybe need for abs())

#### Lagragian unit versus Photo unit

#### Why class based algorithms:

- see https://softwareengineering.stackexchange.com/questions/294076/how-to-represent-an-algorithm-as-a-class