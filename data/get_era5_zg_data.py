# requires configuration file: /home/[yourcomputer]/.cdsapirc
import cdsapi

c = cdsapi.Client()

p_levels = [500]

c.retrieve("reanalysis-era5-pressure-levels", 
    {
    "variable":"z",
    "date":"1979-01-01/2020-03-01",
    "time": "00:00",
    "pressure_level": p_levels,
    "area": [90, 0, 0, 360],
    "product_type": "reanalysis",
    "format": 'netcdf',
    "grid": "1.0/1.0"
    },
"era5_zg_79to21.nc")