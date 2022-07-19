# -*- coding: utf-8 -*-

"""
    Various utility functions needed for the jet-stream metrics and algorithms not belonging to windspeed or spatial utils.
    Built from xarray
"""

# imports
import itertools
import numpy as np
import scipy.signal

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def check_at_least_two_plevs_in_data(data):
    """
    Checks there are at least two pressure-levels (plevs) in xarray dataset

    Parameters
    ----------
    data : xarray.Dataset
        Data to check

    Raises
    ----------
    ValueError :
        If not two pressure levels
    """
    if data["plev"].count() < 2:
        raise ValueError(
            "Need at least 2 pressure levels (plevs) for calculation"
        )


def get_local_minima(arr, axis=0):
    """
    Uses scipy.signal.argrelextrema to get index location of maximum value in array

    Taken from
    https://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array

    Parameters
    ----------
    arr : array-like
        array to find index location of maxima value from
    axis : int
        axis for scipy.signal.argrelextrema
    """
    return scipy.signal.argrelextrema(arr, np.less, axis=axis)


def get_local_maxima(arr, axis=0):
    """
    Uses scipy.signal.argrelextrema to get index location of minimum value in array

    Taken from
    https://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array

    Parameters
    ----------
    arr : array-like
        array to find index location of minima value from
    axis : int
        axis for scipy.signal.argrelextrema
    """
    return scipy.signal.argrelextrema(arr, np.greater, axis=axis)


def is_djf(month):
    """
    Mask used for getting DJF
    """
    return (month == 12) | (month >= 1) & (month <= 2)


def remove_duplicates(arr):
    """
    Removes duplicates from array. From:
    https://stackoverflow.com/questions/2213923/removing-duplicates-from-a-list-of-lists

    Parameters
    ----------
    arr : array-like
        arr to remove duplicates from

    Returns
    ----------
    arr : array-like
        arr with no duplicates
    """

    arr.sort()
    return list(v for v, _ in itertools.groupby(arr))


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


def get_num_of_decimal_places(num):
    """
    Gets number of decimal places in a float

    Parameters
    ----------
    num : float or int
        input number to get decimal places from

    Returns
    ----------
    decimal_places : int
        number of decimal places
    """
    num = "{:f}".format(num).rstrip("0")
    decimal_places = num[::-1].find(".")
    if decimal_places < 0:
        decimal_places = 0
    return decimal_places
