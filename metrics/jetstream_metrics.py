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


def koch_et_al_2006(data, ws_threshold=30):
    """
        TODO: check with chris
        TODO: add equation to this doc

        Returns
        ----------
        weighted_average_ws : DataArray
    """
    # Step 1: get all pressure levels in data as list and make sure hPa TODO: what if mbar? 
    all_plevs_hPa = jetstream_metrics_utils.get_all_plev_hPa(data)
    ## Step 2: calculate sum of weighted windspeed
    print('Step 1: Calculate weighted sum...')
    sum_weighted_ws = 0
    for plev, (i,plev_hPa) in zip(data['plev'], enumerate(all_plevs_hPa)):
        if i != 0:
            plev_hPa = plev_hPa - all_plevs_hPa[i-1]
        sum_weighted_ws += ((data.sel(plev=plev)['ua']** 2 + data.sel(plev=plev)['va']**2)**(1/2)) * plev_hPa
    
    ## Step 3: calculate average weighted
    print('Step 2: Calculate weighted average...')
    weighted_average_ws = sum_weighted_ws * (1/(all_plevs_hPa.max() - all_plevs_hPa.min()))

    ## Step 4: apply threshold
    print('Step 3: Apply windspeed threshold of %s m/s...' % (ws_threshold))
    weighted_average_ws = weighted_average_ws.where(weighted_average_ws >= ws_threshold)
    weighted_average_ws = weighted_average_ws.fillna(0.0)
    return weighted_average_ws


def archer_caldeira_2008(data):
    """
        Will calculate only the mass-weighted wind speed
        Similar to Koch et al. 2006 -> "To overcome this problem, we define jet stream properties via integrated quantities, which are more numerically stable and less grid-dependent than are simple maxima and minima."
    """
    return


def woolings_et_al_2010(data, filter_freq=10):
    """
        Follows an in-text description of 4-steps describing the algorithm mof jet-stream identification from Woolings et al. (2010). 
        Will calculate this metric based on data (regardless of pressure level of time span etc.). 
        TODO: Ask Chris about fourier filtering (step 4)
        TODO: Maybe note about using season or not
    """
    ## Step 1
    dims_for_mean = ['lon', 'plev']
    if 'plev' not in data.dims:
        dims_for_mean = ['lon']
    print('Step 1: calculating long and plev mean...')
    mean_data = data.mean(dims_for_mean)
    ## Step 2
    print('Step 2: Applying %s day lancoz filter...' % (filter_freq))
    lanczos_weights = jetstream_metrics_utils.low_pass_weights(61, 1/filter_freq)
    lanczos_weights_arr = xr.DataArray(lanczos_weights, dims=['window'])
    window_cons = mean_data['ua'].rolling(time=len(lanczos_weights_arr), center=True).construct('window').dot(lanczos_weights_arr)
    ## Step 3
    print('Step 3: Calculating max windspeed and latitude where max windspeed found...')
    max_lat_ws = np.array(list(map(jetstream_metrics_utils.get_latitude_and_speed_where_max_ws, window_cons[:])))
    ## Step 4 — fourier filtering  TODO
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
    

def manney_et_al_2011(data):
    """
        Write function description

        Used in Manney 2011, 2014, 2017 and 2018
    """
    return


def screen_and_simmonds_2013(data):
    """
        Write function description
        Also used in screen and simmonds 2014
        TODO: ask Chris about interpolation method
        TODO: insure that the Earth sphericity is accounted for in the perimeter calculation
    """
    return


def kuang_et_al_2014(data, ws_threshold=30):
    """
        Looks to get event-based jet occurrence percentage and jet center occurrence of (UT)JS
        May take a long time for a lot of data
        TODO: ask chris to check
    """
    print('Step 1. Calculate resultant wind vector')    
    ws_data = jetstream_metrics_utils.get_resultant_wind(data['ua'], data['va'])
    
    print('Step 2. Get values with windspeeds above %s m/s (\'jet-stream occurence\' points)' % (ws_threshold))
    jet_occurence = ws_data.where(ws_data.variable >= ws_threshold)
    
    print('Step 3. Get values of jet-stream centre points (\'jet-stream centre\' points)')
    jet_centres = jetstream_metrics_utils.get_jet_centre_data(jet_occurence)
    return jet_occurence, jet_centres



def francis_vavrus_2015(data, lat_min=20, lat_max=80):
    """
        Write function description
        TODO: maybe add seasonal anomaly calculation 
    """
    ## Step 1 calculate MCI index for data
    print('calculating Meridional Circulation Index from data')
    mci_data = jetstream_metrics_utils.meridional_circulation_index(data) 
    
    ## maybe TODO: Step ?? Calculate anomaly from season
    
    return mci_data



def local_wave_activity(data):
    """
        Introduced by Huang and Nakamura for Potential Vorticity, but then used by:
        Martineau 2017, Chen 2015 and Blackport & Screen 2020 use LWA with 500 hPa zg instead of pv
        TODO: Ask Chris about equation in Blackport 2020 and others
    """
    return


def cattiaux_et_al_2016(data):
    """
        Write function description
        TODO: waiting on Chris' help with the interpolation method from screen and simmonds
    """
    return


def grise_et_al_2017(data):
    """
        Write function description
        See also Ceppi et al. 2012
    """
    return


def ceppi_et_al_2018(data):
    """
        Write function description
        TODO: the centroid??
        TODO: calc for each time period
        "similar methods used in: Chen et al. 2008; Ceppi et al. 2014"
    """
    # lat_mean_vals = []
    # for lat in data['lat']:
    #     lat_mean_vals.append([float(lat), float(data.sel(lat=lat)['ua'].mean()/data['ua'].mean())])
    return


def kern_et_al_2018(data):
    """
        Write function description
        TODO: ask Chris about equation
    """
    return


def simpson_et_al_2018(data):
    """
        Write function description
    """
    return


def lee_et_al_2019(data):
    """
        Write function description
    """
    return


def chemke_and_ming_2020(data):
    """
        Write function description
    """
    return

