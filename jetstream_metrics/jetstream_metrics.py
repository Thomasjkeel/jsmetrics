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


def Koch_et_al_2006(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def woolings_et_al_2010(data, filter_freq=10, lat_min=15, lat_max=75):
    """
        Follows an in-text description of 5-steps describing the algorithm mof jet-stream identification from Woolings et al. (2010). 
        Will calculate this metric based on data (regardless of pressure level of time span etc.). 
        TODO: ask Chris about fourier filtering 
        TODO: maybe note about using season or not
    """
    ## Step 1
    mean_data = data.mean(['lon','plev'])
    ## Step 2
    mean_data = mean_data.sel(lat=slice(lat_min, lat_max))
    ## Step 3
    lanczos_weights = jetstream_metrics_utils.low_pass_weights(61, 1/filter_freq)
    lanczos_weights_arr = xr.DataArray(lanczos_weights, dims=['window'])
    window_cons = mean_data['ua'].rolling(time=len(lanczos_weights_arr), center=True).construct('window').dot(lanczos_weights_arr)
    ## Step 4
    max_lat_ws = np.array(list(map(jetstream_metrics_utils.get_latitude_and_speed_where_max_ws, window_cons[:])))
    # Step 5 — fourier filtering  TODO
    ## make seasonal (DJF, MAM, JJA, SON)
    # seasonal_data = window_cons.resample(time='Q-NOV').mean()
    # filled_seasonal_data = seasonal_data[:,0].fillna(0)
    ## loop over each latitude and calculate high frequency
    # for lat in seasonal_data['lat']:
    #     lat_data = seasonal_data.sel(lat=lat)
    #     lat_data = np.array(lat_data.fillna(0))
    #     filtered_sig = jetstream_metrics_utils.fourier_filter(lat_data)
    #     # code to put back into xarray format
    return max_lat_ws
    

def Hudson_2012(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Screen_and_Simmonds_2014(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Kuang_et_al_2014(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Huang_and_Nakamura_2015(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Cattiaux_et_al_2016(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Ceppi_et_al_2018(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Kern_et_al_2018(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Manney_et_al_2018(data):
    """
        Write function description

        Extend Manney 2011, 2014 and 2017
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Simpson_et_al_2018(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Lee_et_al_2019(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return


def Chemke_and_Ming_2020(data):
    """
        Write function description
    """
    # i.e. will calculate metric based on data (regardless of pressure level of time span etc.)
    return

