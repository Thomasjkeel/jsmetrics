# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to
    identify or classify jet-stream in the literature.

    This file is in order of publish year of the metrics
    (see details_for_all_metrics.py)
"""

# imports
# import collections
import numpy as np

# import matplotlib.pyplot
# import xarray as xr
# import scipy.fftpack
# import scipy.interpolate
# # import shapely.geometry
# from . import data_utils, spatial_utils, windspeed_utils

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def get_sum_weighted_ws(data, all_plevs_hPa):
    """
    Get sum of weighted windspeed.
    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    sum weighted windspeed = integral(p2, p1)(u^2+v^2)^(1/2)dp
    where p1, p2 is min, max pressure level

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind
    all_plevs_hPa : array-like
        list of hPa unit pressure levels

    Returns
    ----------
    sum_weighted_ws : xarray.Dataset
        Data containing sum weighted windspeed values
    """
    if "plev" not in data.coords:
        raise KeyError("Data does not contain coord: 'plev'")

    if not isinstance(all_plevs_hPa, (list, np.ndarray)):
        raise TypeError(
            "array of pressure level needs to be list or numpy.array"
        )

    sum_weighted_ws = 0
    for plev, (i, plev_hPa) in zip(data["plev"], enumerate(all_plevs_hPa)):
        if i != 0:
            plev_hPa = plev_hPa - all_plevs_hPa[i - 1]
        sum_weighted_ws += (
            (data.sel(plev=plev)["ua"] ** 2 + data.sel(plev=plev)["va"] ** 2)
            ** (1 / 2)
        ) * plev_hPa
    return sum_weighted_ws


def get_weighted_average_ws(sum_weighted_ws, all_plevs_hPa):
    """
    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    weighted average windspeed = 1/(p2-p1) * sum average windspeed
    where p1, p2 is min, max pressure level

    Parameters
    ----------
    sum_weighted_ws : xarray.Dataset
        Data containing sum weighted windspeed values
    all_plevs_hPa : array-like
        list of hPa unit pressure levels

    Returns
    ----------
    weighted_average_ws : xarray.Dataset
        Data containing weighted average windspeed values
    """
    if not isinstance(all_plevs_hPa, (list, np.ndarray)):
        raise TypeError(
            "array of pressure level needs to be a list or numpy.array"
        )
    weighted_average_ws = sum_weighted_ws * (
        1 / (all_plevs_hPa.max() - all_plevs_hPa.min())
    )
    return weighted_average_ws


def get_all_hPa_list(data):
    """
    Will get a list of all the pressure levels in the data in hPa/mbar

    Parameters
    ----------
    data : xarray.Dataset
        data with plev coord

    Returns
    ----------
    plev : np.array
        arr containing pressure level list in data in hPa/mbar
    """
    if "plev" not in data.coords:
        raise KeyError("Data does not contain coord: 'plev'")

    if (
        data["plev"].units != "Pa"
        and data["plev"].units != "hPa"
        and data["plev"].units != "mbar"
    ):
        raise ValueError("Plev units need to be mbar, Pa or hPa")

    plevs = np.array([plev for plev in data["plev"]])
    if data["plev"].units == "Pa":
        plevs = plevs / 100
    return plevs
