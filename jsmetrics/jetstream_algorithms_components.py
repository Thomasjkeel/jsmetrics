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
