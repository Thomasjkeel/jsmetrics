import cdsapi

c = cdsapi.Client()

## For 150 - 950 hPa:
# p_levels = [str(z)  for z in ([1] + list(range(150, 1000, 50)))]
p_levels = [700, 850, 925]

c.retrieve("reanalysis-era5-pressure-levels", 
    {
    "variable":"u",
    "date":"1979-01-01/2020-03-01",
    "time": "00:00",
    "pressure_level": p_levels,
    "area": [75, -360, 15, -0],
    "product_type": "reanalysis",
    "format": 'netcdf',
    "grid": "1.0/1.0"
    },
"era5_u_wind_79to21.nc")