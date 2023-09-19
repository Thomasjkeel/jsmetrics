# -*- coding: utf-8 -*-

"""
    Wrappers for checking inputs to the metrics
"""
import functools
import xarray


def sort_xarray_data_coords(coords, ascending=True):
    """
    Sort the coordinates of the input data (e.g. lat or lon). Sort in ascending order by default.
    Used in all jet statistics and jet core algorithms included in this package.

    Parameters
    ----------
    coords : array-like
        Coords to sort in xarray data
    ascending : Boolean (default=True)
        Whether the coords should be sorted in ascending order
    """

    def wrap(func):
        @functools.wraps(func)
        def wrapped_func(*args, **kwargs):
            if "data" in kwargs:
                data = kwargs["data"]
            elif args:
                #  Assumes first argument is data
                data = args[0]
            else:
                return func(*args, **kwargs)
            check_data_is_xarray(data)
            for coord in coords:
                assert (
                    coord in data.coords
                ), f"'{coord}' is not in data. Please check your data and the variable names. Should be like: 'lat', 'lon', 'plev', etc."
                if data[coord].size == 1 and not data[coord].shape == (1,):
                    data = data.expand_dims(
                        coord
                    )  # expand dimensions so the sortby function works if only one in dim
                data = data.sortby(coord, ascending=ascending)
            return func(data, *args[1:], **kwargs)

        return wrapped_func

    return wrap


def check_data_is_xarray(data):
    assert isinstance(data, xarray.DataArray) or isinstance(
        data, xarray.Dataset
    ), "'data' needs to be a xarray DataArray or Dataset type"
