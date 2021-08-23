# -*- coding: utf-8 -*-

"""
    Base classes and function for dealing with windspeed vectors and slices of windspeed
    TODO: Probably needs a better name
"""

### imports
import xarray as xr
import numpy as np

### docs
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


class WindSpeedSlice():
    """
        Base class for windspeed slice
    """
    def __init__(self, data, req_variables=('ua', 'va')):
        ## these will check that the correct data variables and coords are available in the xarray data input
        self.req_variables = req_variables
        self._check_input_data_can_be_used_for_windspeed_slice(data)
        self.values = self._calc_windspeed(data)
        self.values = self.values.rename('ws').to_dataset()
    
    def __init_subclass__(cls, req_coords, *a, **kw):
        cls.req_coords = req_coords

    def __getitem__(self, item):
         return self.values[item]
    
    def __add__(self, other):
        if self.req_coords == other.req_coords and self.req_variables == other.req_variables:
            return xr.concat([self.values, other.values], dim='time')
    
    @staticmethod
    def _calc_windspeed(data):
        """
            The reason this exists is that their may be a point in time where constraitns on this calculation will need to be made i.e. if too much data or too much RAM in use
        """
        return get_resultant_wind(data['ua'], data['va'])
        
    def _check_input_data_can_be_used_for_windspeed_slice(self, data):
        check_if_xarray_dataset_or_array(data)
        check_coord_in_data(data, self.req_coords)
        check_var_in_data(data, self.req_variables)
        check_only_required_coords_are_in_data(data, self.req_coords, to_remove=('bnds',))
    
    def label_slice(self, condition, label):
        return self.values.where(condition, other=label)
    
    def get_values(self):
        return self.values

        
class PressureLevelWindSpeedSlice(WindSpeedSlice, req_coords=('lat', 'lon')):
    """
        Data will be lon*lat
    """ 
    def print_lats(self):
        print(self['lat'])
    
        
class LatitudeWindSpeedSlice(WindSpeedSlice, req_coords=('lat', 'plev')):
    """
        Data will be lon*plev
    """ 
    def print_lats(self):
        print(self['lat'])

        
## Checks for windspeed slices
def check_if_xarray_dataset_or_array(data):
    if not isinstance(data, xr.Dataset) or  isinstance(data, xr.DataArray):
        raise TypeError("input needs to be xarray.DataSet or xarray.DataArray")


def check_coord_in_data(data, req_coords):
    for coord in req_coords:
        if not coord in data.coords:
            raise KeyError("\'%s\' is not in the data" % (coord,))

def check_var_in_data(data, req_variables):
    for var in req_variables:
        if not var in data.variables:
            raise KeyError("\'%s\' is not the data" % (var,))


def check_only_required_coords_are_in_data(data, req_coords, to_remove=()):
    dims = set(data.dims)
    for rem in to_remove:
        try:
            dims.remove(rem)
        except:
            pass
            # print('cannot remove %s from dims' % (rem))

    req_coords = set(req_coords)
    difference = dims.difference(set(req_coords))
    if len(difference) > 0:
        raise ValueError("Unwanted coords in data: %s. Please subset/remove so that the slice can be 2D." % (difference,))
    