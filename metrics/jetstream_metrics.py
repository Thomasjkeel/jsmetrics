# -*- coding: utf-8 -*-

"""
    Metrics used to identify or classify jet-stream in the literature

    All functions should return a xarray.Dataset unless otherwise stated
"""

### imports
import numpy as np
from numpy.core.fromnumeric import mean
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
    print('Step 1: Calculate weighted sum...')
    # Step 1.1: get all pressure levels in data as list and make sure hPa TODO: what if mbar? 
    all_plevs_hPa = jetstream_metrics_utils.get_all_plev_hPa(data)
    # Step 1.2 get weighted sum windspeed
    sum_weighted_ws = jetstream_metrics_utils.get_sum_weighted_ws(data, all_plevs_hPa)
    ## Step 2: calculate average weighted
    print('Step 2: Calculate weighted average...')
    weighted_average_ws = jetstream_metrics_utils.get_weighted_average_ws(sum_weighted_ws, all_plevs_hPa)
    ## Step 3: apply threshold
    print('Step 3: Apply windspeed threshold of %s m/s...' % (ws_threshold))
    weighted_average_ws = weighted_average_ws.where(weighted_average_ws >= ws_threshold)
    weighted_average_ws = weighted_average_ws.fillna(0.0)
    ## turn into dataset
    weighted_average_ws = weighted_average_ws.rename('weighted_average_ws').to_dataset()
    return weighted_average_ws


def archer_caldeira_2008(data):
    """
        Will calculate only the mass-weighted wind speed
        Similar to Koch et al. 2006 -> "To overcome this problem, we define jet stream properties via integrated quantities, which are more numerically stable and less grid-dependent than are simple maxima and minima."
    """
    return data


def woolings_et_al_2010(data, filter_freq=10, window_size=61):
    """
        window_size int: number of days
        Follows an in-text description of 4-steps describing the algorithm mof jet-stream identification from Woolings et al. (2010). 
        Will calculate this metric based on data (regardless of pressure level of time span etc.). 
        TODO: Ask Chris about fourier filtering (step 4)
        TODO: Maybe note about using season or not

        returns:
            max_lat_ws (numpy.ndarray):
    """
    ## Step 1
    print('Step 1: calculating long and/or plev mean...')
    mean_data = jetstream_metrics_utils.get_zonal_mean(data)
    ## Step 2
    print('Step 2: Applying %s day lancoz filter...' % (filter_freq))
    filtered_mean_data = jetstream_metrics_utils.apply_lancoz_filter(mean_data, filter_freq, window_size)
    ## Step 3
    print('Step 3: Calculating max windspeed and latitude where max windspeed found...')
    max_lat_ws = np.array(list(map(jetstream_metrics_utils.get_latitude_and_speed_where_max_ws, filtered_mean_data[:])))
    ## Step 4 — fourier filtering  TODO
    ### make seasonal (DJF, MAM, JJA, SON)
    # seasonal_data = filtered_mean_data.resample(time='Q-NOV').mean()
    # filled_seasonal_data = seasonal_data[:,0].fillna(0)
    ## loop over each latitude and calculate high frequency
    # for lat in seasonal_data['lat']:
    #     lat_data = seasonal_data.sel(lat=lat)
    #     lat_data = np.array(lat_data.fillna(0))
    #     filtered_sig = jetstream_metrics_utils.fourier_filter(lat_data)
    #     ### code to put back into xarray format

    return max_lat_ws
    
        
def manney_et_al_2011(data, ws_core_threshold=40, ws_boundary_threshold=30):
    """
        Write function description

        Used in Manney 2011, 2014, 2017 and 2018
    """
        
    if data['time'].count() > 1:
        for time_coord in data['time']:
            sub_data = data.sel(time=time_coord)
            core_alg = jetstream_metrics_utils.JetStreamCoreIdentificationAlgorithm(sub_data, ws_core_threshold=ws_core_threshold, ws_boundary_threshold=ws_boundary_threshold)
            yield core_alg
    else:
        core_alg = jetstream_metrics_utils.JetStreamCoreIdentificationAlgorithm(data, ws_core_threshold=ws_core_threshold, ws_boundary_threshold=ws_boundary_threshold)
        yield core_alg
        
        
def screen_and_simmonds_2013(data):
    """
        Write function description
        Also used in screen and simmonds 2014
        TODO: ask Chris about interpolation method
        TODO: insure that the Earth sphericity is accounted for in the perimeter calculation
    """
    return


def kuang_et_al_2014(data, occurence_ws_threshold=30):
    """
        Looks to get event-based jet occurrence percentage and jet center occurrence of (UT)JS
        May take a long time for a lot of data
        TODO: ask chris to check
    """ 
    if data['time'].count() > 1:
        for time_coord in data['time']:
            sub_data = data.sel(time=time_coord)
            occurence_alg = jetstream_metrics_utils.JetStreamOccurenceAndCentreAlgorithm(sub_data, occurence_ws_threshold=occurence_ws_threshold)
            yield occurence_alg
    else:
        occurence_alg = jetstream_metrics_utils.JetStreamOccurenceAndCentreAlgorithm(data, occurence_ws_threshold=occurence_ws_threshold)
        yield occurence_alg



def francis_vavrus_2015(data, lat_min=20, lat_max=80):
    """
        Write function description
        TODO: maybe add seasonal anomaly calculation 
    """
    ## Step 1 calculate MCI index for data
    print('calculating Meridional Circulation Index from data')
    mci_data = jetstream_metrics_utils.meridional_circulation_index(data) 
    
    ## maybe TODO: Step ?? Calculate anomaly from season
    mci_data = mci_data.rename('mci').to_dataset()

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
        TODO: work out if relevant as this method also uses poleward edge of sub-tropical dry zone and poleward edge of Hadley cell derived from precip. record 
        TODO: and add to dict if relevant
    """
    return


def ceppi_et_al_2018(data):
    """
        Write function description
        TODO: what is meant by the centroid??
        "similar methods used in: Chen et al. 2008; Ceppi et al. 2014"

        Returns: centroid latitude of u-wind for one day
    """
    if data['time'].count() > 1:
        for time_coord in data['time']:
            sub_data = data.sel(time=time_coord)
            yield jetstream_metrics_utils.get_centroid_jet_lat(sub_data)
    else:
        yield jetstream_metrics_utils.get_centroid_jet_lat(data)


def kern_et_al_2018(data):
    """
        Write function description
        TODO: ask about equation
    """
    return


def simpson_et_al_2018(data):
    """
        Write function description
        TODO: ask about interpolation
        Before comparing the variability between the reanalyses and the models, each dataset is first interpolated onto a 2*2
        longitude–latitude grid using a cubic spline interpolation and then isotropically smoothed in the spectral domain retaining
         only scales larger than total wavenumber 42 according to Sardeshmukh and Hoskins [1984, their Eq. (9) with n0=42 and r=1]. T
    """
    return


def lee_et_al_2019(data):
    """
        Write function description
        TODO: work out if relevant as wind shear of jet-stream rather than location
        variable used: temperature, u-wind, specific gas constant for dry air, the Coriolis parameter, pressure, northward distance.
        TODO: and add to dict if relevant
    """
    return


def chemke_and_ming_2020(data):
    """
        Write function description
        TODO: ask about equation
    """
    return

