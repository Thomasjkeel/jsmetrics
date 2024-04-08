# -*- coding: utf-8 -*-

"""
    Various utility functions needed for the jet-stream metrics and algorithms not belonging to windspeed or spatial utils.
    The module is built from xarray data-structures, so this file contains all the stuff that helps the library handle xarray dataset/dataarrays

    Classes and Functions ordered alphabetically
"""

# imports
import cftime
import itertools
import numpy as np
import xarray as xr
import scipy.signal

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def add_num_of_days_to_360Datetime(datetime_360day, num_of_days_to_add):
    """
    Adds a number of days to cftime.Datetime360Day format

    Parameters
    ----------
    datetime_360day : cftime.Datetime360Day
        Date to add days to

    num_of_days_to_add : int or float
        Number of days to add to 360 day format

    Returns
    ----------
    new_360day_date : cftime.Datetime360Day
        Date with added number of days

    Raises
    ----------
    AssertionError :
        If not '360_day' format or not a cftime.Datetime object
    ValueError :
        If number of days not more than 0
    """
    assert hasattr(
        datetime_360day, "calendar"
    ), "date type inputted is not in cftime format."
    assert (
        getattr(datetime_360day, "calendar") == "360_day"
    ), "input date is not '360_day' format."
    if num_of_days_to_add <= 0:
        raise ValueError("Number of days to add to date is not more than 0")

    new_day = ((datetime_360day.day + num_of_days_to_add) - 30) % 30
    if new_day == 0:
        new_day = 30
    num_of_months_to_add = (datetime_360day.day + num_of_days_to_add) / 30
    if num_of_months_to_add % 1 == 0:
        num_of_months_to_add = num_of_months_to_add - 1
    else:
        num_of_months_to_add = np.floor(num_of_months_to_add)
    new_month = ((datetime_360day.month + num_of_months_to_add) - 12) % 12
    if new_month == 0:
        new_month = 12
    num_of_years_to_add = (datetime_360day.month + num_of_months_to_add) / 12
    if num_of_years_to_add % 1 == 0:
        num_of_years_to_add = num_of_years_to_add - 1
    else:
        num_of_years_to_add = np.floor(num_of_years_to_add)

    new_year = datetime_360day.year + num_of_years_to_add
    new_360day_date = cftime.Datetime360Day(
        day=new_day,
        month=new_month,
        year=new_year,
        hour=datetime_360day.hour,
        minute=datetime_360day.minute,
        second=datetime_360day.second,
    )
    return new_360day_date


def add_num_of_days_to_NoLeapDatetime(datetime_noleap, num_of_days_to_add):
    """
    Adds a number of days to cftime.DatetimeNoLeap format

    Parameters
    ----------
    datetime_noleap : cftime.DatetimeNoLeap
        Date to add days to

    num_of_days_to_add : int or float
        Number of days to add to noleap day format

    Returns
    ----------
    new_noleap_date : cftime.DatetimeNoLeap
        Date with added number of days

    Raises
    ----------
    AssertionError :
        If not 'noleap' format or not a cftime.Datetime object
    ValueError :
        If number of days not more than 0
    """
    assert hasattr(
        datetime_noleap, "calendar"
    ), "date type inputted is not in cftime format."
    assert (
        getattr(datetime_noleap, "calendar") == "noleap"
    ), "input date is not 'no_leap' format."
    if num_of_days_to_add <= 0:
        raise ValueError("Number of days to add to date is not more than 0")

    return cftime.DatetimeNoLeap.fromordinal(
        cftime.DatetimeNoLeap.toordinal(datetime_noleap) + num_of_days_to_add,
        calendar="noleap",
        has_year_zero=True,
    )


def add_pad_to_array(arr, pad_width=1, constant_values=0):
    """
    Add a edge of constant values around a numpy array

    Parameters
    ----------
    arr : array-like
        A 2-D array with numeric values (i.e. dtypes float or int)
    pad_width : int
        Width of values to add to edges (default: 1)
    constant_values : int
        Value to add to edges (default: 0)

    Returns
    ----------
    padded_arr : array-like
        A 2-D array with new dimensions as a pad has been added at the edges
    """
    padded_arr = np.pad(
        arr, pad_width, mode="constant", constant_values=constant_values
    )
    return padded_arr


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


def check_plev_units(data, expected_plev_units):
    """
    Will check that pressure level units are in data plev coord (i.e. Pa, hPa, mbar, millibars).

    Parameters
    ----------
    data : xr.DataArray or xr.Dataset
        Data containing plev coord and units for plev
    expected_plev_units : list
        List of names of plev units (i.e. mbar, Pa, hPa, etc.)

    Returns
    ----------
    units : str
        Units of plev coord

    Raises
    ----------
    KeyError
        If plev is not a coord in input 'data'. And if the units of plev are not in expected plev units

    """
    if "plev" not in data.coords:
        raise KeyError(
            "Please provide a plev coordinate for data to determine plev units"
        )
    if not hasattr(data["plev"], "units"):
        raise KeyError(
            f'You will need assign units (i.e {expected_plev_units}) to plev to run this method e.g. data["plev"] = data["plev"].assign_attrs(units="Pa")'
        )

    if data["plev"].units not in expected_plev_units:
        raise KeyError(
            f'plev unit: {data["plev"].units} is not in {expected_plev_units}. You will need assign units to plev to run this method  e.g. data["plev"] = data["plev"].assign_attrs(units="Pa")'
        )
    else:
        return data["plev"].units


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
    idx = np.argmin(np.abs(array - value))
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


def filter_local_extremes_to_min_distance(local_extrema, min_distance_threshold=2):
    """
    Filter local extremes (i.e. outputs of minima or maxima from scipy.signal.argrelextrema)
    so that no neighbours remain in array.

    Parameters
    ----------
    local_extrema : array-like
        likely (stacked) outputs of scipy.signal.argrelextrema
    min_distance_threshold : int
        Minimum distances between indexes in array (default: 2 indexes)

    Returns
    -------
    filtered_extrema : array-like
        (stacked) outputs of scipy.signal.argrelextrema with values less than a min_distance_threshold away removed

    Examples
    --------
    maxima_indices = np.column_stack(jsmetrics.utils.data_utils.get_local_maxima(current['ws'].data))
    filtered_indices = filter_local_extremes_so_no_neighbours(maxima_indices, min_distance_threshold=2)
    """
    filtered_extrema = []

    for idx in local_extrema:
        if not any(
            np.linalg.norm(idx - existing_maxima) < min_distance_threshold
            for existing_maxima in filtered_extrema
        ):
            filtered_extrema.append(idx)

    filtered_extrema = np.array(filtered_extrema)
    return filtered_extrema


def find_intersection_between_two_array_of_arrays(array1, array2):
    """
    Find the intersection between two arrays of arrays. See examples for example.

    Parameters
    ----------
    array1 : np.array
        First array of arrays to compare with array2
    array2 : np.array
        Second array of arrays to compare with array1

    Returns
    -------
    intersection : np.array
        Intesection of the two arrays returning only the arrays within the original two that are in both.

    Examples
    --------
    array1 = [[1, 2], [3, 4]]
    array2 = [[1, 2], [4, 3]]
    find_intersection_between_two_array_of_arrays(array1, array2)
    # returns [[1, 2]]

    """
    intersection = []
    for arr1 in array1:
        for arr2 in array2:
            if np.array_equal(arr1, arr2):
                intersection.append(arr1)
                break
    intersection = np.array(intersection)
    return intersection


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


def remove_unwanted_coords_from_data(
    data, wanted_coords, unwanted_coords=(), show_error=False
):
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
            if show_error:
                print(e)

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
    Rescale latitude resolution to input lat resolution in degrees.

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani 2014 https://doi.org/10.1002/2013GL058466
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


def slice_array_by_index_breaks(array_to_slice, index_breaks):
    """
    Will break an array down into segments based on a list of indexes containing information about where to create slices

    Parameters
    ----------
    array_to_slice : array-like
        Array to break down based on slice breaks
    index_breaks : array-like
        Indexes from which to slice array

    Returns
    ----------
    output: list
        Broken down slices of original array_to_slice
    """
    return [
        array_to_slice[i:j] for i, j in zip([0] + index_breaks, index_breaks + [None])
    ]
