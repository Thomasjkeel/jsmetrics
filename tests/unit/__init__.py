import numpy as np
import pandas as pd
import xarray as xr
import random


def set_up_test_uv_data():
    u_data = xr.open_dataset(
        "tests/data/ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20150105.nc"
    )
    v_data = xr.open_dataset(
        "tests/data/va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20150105.nc"
    )
    data = xr.merge([u_data, v_data])
    return data


def set_up_test_u_data():
    data = xr.open_dataset(
        "tests/data/ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20150105.nc"
    )
    return data


def set_up_test_zg_data():
    data = xr.open_dataset(
        "tests/data/zg_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20150105.nc"
    )
    return data


def make_fake_data(data, variable_name):
    random.seed(42)
    data[variable_name] = data[variable_name] + random.randint(0, 5)
    return data


def make_fake_seasonal_data(data):
    data["time"] = np.array(
        [
            "2015-01-01T00:00:00.000000000",
            "2015-02-01T00:00:00.000000000",
            "2015-03-01T00:00:00.000000000",
            "2015-04-01T00:00:00.000000000",
            "2015-05-01T00:00:00.000000000",
        ],
        dtype="datetime64[ns]",
    )
    return data


def set_up_nan_dataset():
    lon = [[99.32, 99.83], [99.23, 99.73]]
    lat = [[42.25, 42.21], [42.63, 42.59]]
    time = pd.date_range("2014-09-06", periods=3)
    reference_time = pd.Timestamp("2014-09-05")
    an_array = np.empty((2, 2, 3))
    an_array[:] = np.nan
    da = xr.DataArray(
        data=an_array,
        dims=["x", "y", "time"],
        coords=dict(
            lon=(["x", "y"], lon),
            lat=(["x", "y"], lat),
            time=time,
            reference_time=reference_time,
        ),
    )
    return da
