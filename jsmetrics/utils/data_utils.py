# -*- coding: utf-8 -*-

"""
    Various utility functions needed for the jet-stream metrics and algorithms not belonging to windspeed or spatial utils.
    The module is built from xarray data-structures, so this file contains all the stuff that helps the library handle xarray dataset/dataarrays

    Classes and Functions ordered alphabetically
"""

# imports
import itertools
import numpy as np
import xarray as xr
import scipy.signal

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def check_at_least_n_plevs_in_data(data, n_plevs):
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
    if data["plev"].count() < n_plevs:
        raise ValueError(
            "Need at least %s pressure levels (plevs) for calculation" % (n_plevs)
        )


def check_coords_in_data(data, req_coords):
    """
    Will check if a given coordinates are in the data and raises a key error if any of them are not.
    This function is needed to check before algorithm continues and uses too much memory.
    Built from xarray.

    Parameters
    ----------
    data : xarray.Dataset
        Data to check
    req_coords : array-like or tuple
        Coordinates to check if in data

    Raises
    ----------
    KeyError :
        If any given coord is not in data
    """
    for coord in req_coords:
        if coord not in data.coords:
            raise KeyError("'%s' is not in the data" % (coord,))


def check_if_data_is_xarray_datatype(data):
    """
    What it says on the tin.

    Parameters
    ----------
    data : xarray.Dataset
        Data to check

    Raises
    ----------
    TypeError :
        If input is not an xarray data type
    """
    if not isinstance(data, xr.Dataset) or isinstance(data, xr.DataArray):
        raise TypeError("input needs to be xarray.DataSet or xarray.DataArray")


def check_variables_in_data(data, req_variables):
    """
    What it says on the tin.
    Built from xarray

    Parameters
    ----------
    data : xarray.Dataset
        Data to check
    req_variables : array-like
        Variables needed in data

    Raises
    ----------
    KeyError :
        If variables not in data
    """
    for var in req_variables:
        if var not in data.variables:
            raise KeyError("'%s' is not the data" % (var,))


def find_nearest_value_to_array(array, value):
    """
    Will find the nearest value to a given array
    Built for use with xarray
    adapted from: https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array

    Parameters
    ----------
    array : array-like
        Array to reference
    value : numeric-like
        Number to check nearest value in input array

    Returns
    ----------
    output : numeric
        Closest value to input value in input array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


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
    return decimal_places


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


def remove_unwanted_coords_from_data(data, wanted_coords, unwanted_coords=()):
    """
    What it says on the tin.
    Built from xarray

    Parameters
    ----------
    data : xarray.Dataset
        Data to check

    wanted_coords : array-like or tuple
        Coords to retain in data
    unwanted_coords : array-like or tuple
        Coords to remove from data


    Raises
    ----------
    ValueError :
        If coord cannot be removed from data
    """
    dims = set(data.dims)
    for rem in unwanted_coords:
        try:
            dims.remove(rem)
        except Exception as e:
            # a little sloppy, but probably okay
            e

    wanted_coords = set(wanted_coords)
    difference = dims.difference(set(wanted_coords))
    if len(difference) > 0:
        raise ValueError(
            "Unwanted coords in data: %s.\
             Please subset/remove so that the slice can be 2D."
            % (difference,)
        )


def rescale_lat_resolution(lats, lat_resolution):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani 2017 https://doi.org/10.1175/JCLI-D-16-0849.1
    & Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1


    Parameters
    ----------
    lats : xr.DataArray or array-like
        Array of latitude values
    lat_resolution : int or float
        Latitude resolution in degrees

    Returns
    ----------
    output : numpy.array
        Rescaled array of latitude values
    """
    return np.arange(min(lats), max(lats) + lat_resolution, lat_resolution)
