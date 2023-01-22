# -*- coding: utf-8 -*-

"""
    Operations needed for the jet-stream metrics and jet-stream algorithms that specifically operate on windspeed data.
    Includes the base classes and function for dealing with windspeed vectors and lat/lon or lat/plev slices of windspeed (so called: WindSpeedSlice class)

    Classes and Functions ordered alphabetically.
"""

# imports
import numpy as np
from . import data_utils

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def get_resultant_wind(u, v):
    """
    Gets wind vector from u-wind and v-wind
    """
    return np.sqrt(u**2 + v**2)


def get_wind_direction_in_degrees(u, v):
    """
    Gets wind direction from u-wind and v-wind
    In degrees (0-360)
    """
    return (180 + (180 / np.pi) * np.arctan2(u, v)) % 360


def get_zonal_mean(data):
    """
    Will get the zonal mean either by pressure level (plev) or for one layer

    Parameters
    ----------
    data : xarray.Dataset
        Data containing lon and plev coords

    Returns
    ----------
    zonal_mean : xarray.DataSet
        zonal mean data

    Raises
    ----------
    KeyError
        when 'lon' not discovered as coord
    """
    if "lon" not in data.coords:
        raise KeyError("data does not contain 'lon' coord")

    coords_for_mean = ["lon", "plev"]
    if "plev" not in data.coords or int(data["plev"].count()) == 1:
        coords_for_mean = ["lon"]
    zonal_mean = data.mean(coords_for_mean)
    return zonal_mean


class WindSpeedSlice:
    """
    Base class for windspeed slice.

    Methods
    -------
    label_slice:
        Label the windspeed slice using where condtion
    get_values:
        Get values using xarray
    """

    def __init__(self, data, req_variables=("ua", "va")):
        # these will check that the correct data variables/coords
        # are available in the xarray data input
        self.req_variables = req_variables
        self._check_input_data_can_be_used_for_windspeed_slice(data)
        self.values = self._calc_windspeed(data)
        self.values = self.values.rename("ws").to_dataset()

    def __init_subclass__(cls, req_coords, *a, **kw):
        cls.req_coords = req_coords

    def __getitem__(self, item):
        return self.values[item]

    @staticmethod
    def _calc_windspeed(data):
        """
        The reason this exists is that their may be a point in time where
        constraitns on this calculation will need to be made
        i.e. if too much data or too much RAM in use
        """
        return get_resultant_wind(data["ua"], data["va"])

    def _check_input_data_can_be_used_for_windspeed_slice(self, data):
        """
        Check input data is valid for the windspeed slice.
        """
        data_utils.check_if_data_is_xarray_datatype(data)
        data_utils.check_coords_in_data(data, self.req_coords)
        data_utils.check_variables_in_data(data, self.req_variables)
        data_utils.remove_unwanted_coords_from_data(
            data, wanted_coords=self.req_coords, unwanted_coords=("bnds",)
        )

    def label_slice(self, condition, label):
        """
        Returns labels where condition is met.
        """
        return self.values.where(condition, other=label)

    def get_values(self):
        """
        Get values.
        """
        return self.values


class PressureLevelWindSpeedSlice(WindSpeedSlice, req_coords=("lat", "lon")):
    """
    Data will be lon*lat

    Methods
    -------
    label_slice:
        Label the windspeed slice using where condtion
    get_values:
        Get values using xarray
    """


class LatitudeWindSpeedSlice(WindSpeedSlice, req_coords=("lat", "plev")):
    """
    Data will be lon*plev

    Methods
    -------
    label_slice:
        Label the windspeed slice using where condtion
    get_values:
        Get values using xarray
    """
