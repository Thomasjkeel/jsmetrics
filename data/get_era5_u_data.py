# requires configuration file: /home/[yourcomputer]/.cdsapirc
import cdsapi

c = cdsapi.Client()

p_levels = [250, 300, 500, 700, 850, 925]

c.retrieve("reanalysis-era5-pressure-levels", 
    {
    "variable":"u",
    "date":"1981-01-01/2021-07-01",
    "time": "00:00",
    "pressure_level": p_levels,
    "area": [90, 0, 0, 360],
    "product_type": "reanalysis",
    "format": 'netcdf',
    "grid": "1.0/1.0"
    },
"era5_u_wind_81to21.nc")