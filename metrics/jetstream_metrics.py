# -*- coding: utf-8 -*-

"""
    Metrics used to identify or classify jet-stream in the literature
"""

### imports
import numpy as np
import xarray as xr
from .import jetstream_metrics_utils
# import jetstream_metrics.jetstream_metrics_utils

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def koch_et_al_2006(data):
    """
        TODO: Ask Chris about the equation
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def woolings_et_al_2010(data, filter_freq=10, lat_min=15, lat_max=75):
    """
        Follows an in-text description of 5-steps describing the algorithm mof jet-stream identification from Woolings et al. (2010). 
        Will calculate this metric based on data (regardless of pressure level of time span etc.). 
        TODO: Ask Chris about fourier filtering 
        TODO: Maybe note about using season or not
    """
    ## Step 1
    dims_for_mean = ['lon', 'plev']
    if 'plev' not in data.dims:
        dims_for_mean = ['lon']
    print('Step 1: calculating long and plev mean...')
    mean_data = data.mean(dims_for_mean)
    ## Step 2
    print('Step 2: Subsetting to lat to between %s and %s...' % (lat_min, lat_max))
    ### check lat max and min are the correct way around
    if mean_data.lat[0] > mean_data.lat[-1]:
        mean_data = mean_data.reindex(lat=list(reversed(mean_data.lat)))
    mean_data = mean_data.sel(lat=slice(lat_min, lat_max))
    ## Step 3
    print('Step 3: Applying %s day lancoz filter...' % (filter_freq))
    lanczos_weights = jetstream_metrics_utils.low_pass_weights(61, 1/filter_freq)
    lanczos_weights_arr = xr.DataArray(lanczos_weights, dims=['window'])
    window_cons = mean_data['ua'].rolling(time=len(lanczos_weights_arr), center=True).construct('window').dot(lanczos_weights_arr)
    ## Step 4
    print('Step 4: Calculating max windspeed and latitude where max windspeed found...')
    max_lat_ws = np.array(list(map(jetstream_metrics_utils.get_latitude_and_speed_where_max_ws, window_cons[:])))
    ## Step 5 — fourier filtering  TODO
    ### make seasonal (DJF, MAM, JJA, SON)
    # seasonal_data = window_cons.resample(time='Q-NOV').mean()
    # filled_seasonal_data = seasonal_data[:,0].fillna(0)
    ## loop over each latitude and calculate high frequency
    # for lat in seasonal_data['lat']:
    #     lat_data = seasonal_data.sel(lat=lat)
    #     lat_data = np.array(lat_data.fillna(0))
    #     filtered_sig = jetstream_metrics_utils.fourier_filter(lat_data)
    #     ### code to put back into xarray format
    return max_lat_ws
    


def screen_and_simmonds_2014(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def kuang_et_al_2014(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def huang_and_nakamura_2015(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def cattiaux_et_al_2016(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def ceppi_et_al_2018(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def kern_et_al_2018(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def manney_et_al_2018(data):
    """
        Write function description

        Extend Manney 2011, 2014 and 2017
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def simpson_et_al_2018(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def lee_et_al_2019(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def chemke_and_ming_2020(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return

