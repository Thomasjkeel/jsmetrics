# -*- coding: utf-8 -*-

"""
    Statistics and algorithms for determining the "waviness" of upper-level mean flow within a given
    time window. These metrics only have meaning at an integrated global scale.

    The following metrics each return a xarray.Dataset and are ordered by paper publish year.
"""

# imports
import xarray
from . import waviness_metrics_components
from jsmetrics.utils import spatial_utils
from jsmetrics.core.check_data import sort_xarray_data_coords

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


@sort_xarray_data_coords(coords=["lat", "lon"])
def francis_vavrus_2015(data):
    r"""
    This method calculates a waviness metric: Meridional Circulation Index (MCI) from u- and v-components of wind.

    MCI is calculated as follows:

    .. math::
        MCI = \frac{v*|v|}{u^2+v^2}

    When MCI = 0, the wind is purely zonal, and when MCI= 1 (-1), the flow is from the South (North).

    This method was originally introduce in Francis & Vavrus (2015) https://doi.org/10.1088/1748-9326/10/1/014005
    and is described in the Section entitled: 'Meridional circulation Index (MCI)'.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing the one output: 'mci'

    Notes
    -----
    In the original methodology, MCI is expressed in terms of seasonal anomaly, we show you how to do this in 'Examples'

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        uv_data = xr.open_dataset('path_to_uv_data')

        # Subset dataset to range used in original methodology (500 hPa & 20-80 N)):
        uv_sub = uv_data.sel(plev=500, lat=slice(20, 80))

        # Run statistic:
        fv15 = jsmetrics.waviness_metrics.francis_vavrus_2015(uv_sub)

        # Express MCI as a seasonal anomaly
        fv15_seasonal_anomalies = (mci['mci'] - mci['mci'].groupby('time.season').mean())
    """
    #  Step 1. calculating Meridional Circulation Index from data
    data["mci"] = waviness_metrics_components.calc_meridional_circulation_index(data)
    return data


@sort_xarray_data_coords(coords=["lat", "lon"])
def cattiaux_et_al_2016(data):
    r"""
    This method calculates a sinousity metric for upper-air flow using geopotential height.
    The value of sinuosity is selected using an isohypse which precisely corresponds to the Z500 average over 30-70∘N .
    Then this value is compared to the perimeter at 50∘N to calculate sinuosity.

    This method was first introduce in Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309 and
    is described in section 3.1 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'zg', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing the two outputs: 'sinousity' and 'zonal_mean_zg_30Nto70N'

    Notes
    -----
    The original implementation used the R package 'geosphere' to calculate sinuosity. Here we use the Python package
    Shapely to calculate great circle distances between the two perimeters to calculate sinuosity.

    **Moderately slow method:** currently takes a moderate amount of time i.e. 2 seconds per 100 time units\
    with 1 plev on AMD Ryzen 5 3600 6-core processor

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        zg_data = xr.open_dataset('path_to_zg_data')

        # Subset dataset to range used in original methodology (500 hPa & 00-90 N)):
        zg_sub = zg_data.sel(plev=500, lat=slice(0, 90))

        # Run statistic:
        c16 = jsmetrics.waviness_metrics.cattiaux_et_al_2016(zg_sub)
    """
    # Check input is xarray.Dataset
    if isinstance(data, xarray.DataArray):
        data = data.to_dataset()

    #  Step 1. get zonal average for each timestep
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        data = data.expand_dims("time")
    data["zonal_mean_zg_30Nto70N"] = (
        data["zg"].sel(lat=slice(30, 70)).groupby("time").mean(...)
    )

    #  Step 2. Get latitude circle of 50 N
    circle_50N = spatial_utils.get_latitude_circle_linestring(50, 0, 360)

    #  Step 3. Loop over each time step and calculate sinuosity
    if data["time"].size == 1:
        # shrink dims again
        data = data.isel(time=0)
        output = waviness_metrics_components.get_sinuosity_of_zonal_mean_zg(
            data, circle_50N
        )
    else:
        output = data.groupby("time").map(
            lambda row: waviness_metrics_components.get_sinuosity_of_zonal_mean_zg(
                row, circle_50N
            )
        )
    return output
