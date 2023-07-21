# -*- coding: utf-8 -*-

"""
    Jet-stream waviness metrics used in the literature.

    Classes and Functions ordered by paper publish year.
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


# @sort_xarray_data_coords(coords=["lat", "lon"])
# def screen_and_simmonds_2013(data):
#     """
#     Slightly adjusted in Screen and Simmonds 2014
#     Method from Screen & Simmonds (2013) https://doi.org/10.1002/grl.50174

#     Parameters
#     ----------
#     data : xarray.Dataset
#         Data containing geopotential height (zg)

#     Returns
#     ----------
#     data : xarray.Dataset
#         Data containing u- and v- component wind speed
#     """
#     if isinstance(data, xarray.DataArray):
#         data = data.to_dataset()
#     return data


# @sort_xarray_data_coords(coords=["lat", "lon"])
# def local_wave_activity(data):
#     """
#     Introduced by Huang and Nakamura for Potential Vorticity, but then used by:
#     Martineau 2017, Chen 2015 and Blackport & Screen 2020 use LWA
#     with 500 hPa zg instead of pv

#     Parameters
#     ----------
#     data : xarray.Dataset
#         Data containing geopotential height (zg)
#     """
#     return data


@sort_xarray_data_coords(coords=["lat", "lon"])
def francis_vavrus_2015(data):
    """
    Calculates the Meridional Circulation Index (MCI). When MCI = 0, the wind is purely zonal, and when MCI= 1 (-1), the flow is
    from the South (North).

    Method from Francis & Vavrus (2015) https://doi.org/10.1088/1748-9326/10/1/014005

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind

    Returns
    ----------
    data : xarray.Dataset
        Data containing MCI
    """
    #  Step 1. calculating Meridional Circulation Index from data
    data["mci"] = waviness_metrics_components.calc_meridional_circulation_index(data)

    #  Step 2. TODO Calculate anomaly from season
    # maybe TODO: Step 3 Calculate anomaly from season
    return data


@sort_xarray_data_coords(coords=["lat", "lon"])
def cattiaux_et_al_2016(data):
    """
    A sinousity metric for upper-air flow.
    Calculates, for each time unit, the value of the selected isohypse precisely corresponds to the Z500 average over 30–70∘N .
    Then uses the perimeter of this isohype and around 50 .N to calculate sinuosity
    Method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    NOTE: Currently takes a moderate amount of time i.e. 2 seconds per 100 time unit with 1 plev on AMD Ryzen 5 3600 6-core processor

    Parameters
    ----------
    data : xarray.Dataset
        Data containing geopotential height (zg)

    Returns
    ----------
    data : xarray.Dataset
        Data containing sinuosity of zonal mean by each time unit
    """
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
