# -*- coding: utf-8 -*-

"""
    Metrics used to identify or classify jet-stream in the literature
"""

### imports
import numpy as np
import xarray as xr
### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def woolings_et_al_2010(data, filter_freq=10, lat_min=15, lat_max=75):
    """
        Write function description
        TODO: maybe default values for freq and lat/lon bands
        TODO: add a useful exception capture
        TODO: maybe note about using season or not
    """
    # i.e. will calculate woolings metric based on data (regardless of pressure level of time span etc.)
    mean_data = data.mean(['lon','plev'])
    mean_data = mean_data.sel(lat=slice(lat_min, lat_max))
    lanczos_weights = low_pass_weights(61, 1/10)
    lanczos_weights_arr = xr.DataArray(lanczos_weights, dims=['window'])
    window_cons = mean_data['ua'].rolling(time=len(lanczos_weights_arr), center=True).construct('window').dot(lanczos_weights_arr)
    max_lat_ws = np.array(list(map(get_latitude_and_speed_where_max_ws, window_cons[:])))
    return max_lat_ws


def get_latitude_and_speed_where_max_ws(data_row, latitude_col='lat'):
    if not data_row.isnull().all():
        max_speed_loc = data_row.argmax().data
        max_speed = data_row[max_speed_loc]
        lat_at_max = float(max_speed['lat'].values)
        speed_at_max = float(max_speed.data)
        return lat_at_max, speed_at_max 
    else:
        return None, None
    



def Koch_et_al_2006(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Ceppi_et_al_2015(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.
    
    A low-pass filter removes short-term random fluctations in a time series

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    TAKEN FROM: https://scitools.org.uk/iris/docs/v1.2/examples/graphics/SOI_filtering.html

    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[0+(window%2):-1] # edited from w[1:-1]