# -*- coding: utf-8 -*-

"""
    Plotting classes and functions
"""

### imports
import xarray as xr

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


class Plotter():
    """
        Base class for plotter. Input has to be xr.DataArray
    """
    def __init__(self, data):
        assert type(data) == xr.DataArray, "input data to Plotter() class needs be xarray.DataArray type"
        self.data = data
        self._check_required_coords_are_in_data()
        
    def __init_subclass__(cls, req_coords, *a, **kw):
        cls.req_coords = req_coords
        
    def _check_required_coords_are_in_data(self):
        if not hasattr(self.data, 'req_coords'):
            print('Using base Plotter() class')
            self.req_coords = ()
        for coord in self.req_coords:
            assert coord in self.data.dims, "'%s' not in data and is required for Plotter" % (coord)
        

class TimePlotter(Plotter, req_coords=('time',)):
    pass


class TimeLatLonPlotter(Plotter, req_coords=('lat', 'lon')):
    pass
        

    
class TimePlevLonPlotter(Plotter, req_coords=('plev','lon')):
    pass
    