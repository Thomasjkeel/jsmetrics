# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to identify or classify jet-stream in the literature
"""

### imports
import numpy as np
import xarray as xr
from scipy import fftpack
import collections

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def get_resultant_wind(u_data, v_data):
    """
        TODO: work out where and how this will be used
    """
    return np.sqrt( u_data**2 + v_data**2 )



def get_all_plev_hPa(data):
    """
        Will get a list of all the pressure levels in the data in hPa 
    """
    plevs = np.array([plev for plev in data['plev']])
    if data['plev'].units == 'Pa':
        plevs = plevs/100 
    return plevs


def get_latitude_and_speed_where_max_ws(data_row, latitude_col='lat'):
    """
        Write function description
    """
    if not data_row.isnull().all():
        max_speed_loc = data_row.argmax().data
        max_speed = data_row[max_speed_loc]
        lat_at_max = float(max_speed['lat'].values)
        speed_at_max = float(max_speed.data)
        return lat_at_max, speed_at_max 
    else:
        return None, None


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


def fourier_filter(data, timestep=1):
    """
        Carries out a Fourier transform for high frequency filtering
        TAKEN FROM: https://scipy-lectures.org/intro/scipy/auto_examples/plot_fftpack.html
        NOTE: NOT CURRENTLY WORKING PROPERLY
        
        Parameters
        ----------
        data : np.array (1-d) 
            time series data
        timestep : float or int
            number used in the Discrete Fourier Transform sample frequencies (fftfreq)
    """
    fourier_transform = fftpack.fft(data)
    
    # The corresponding frequencies TODO: what does this do?
    sample_freq = fftpack.fftfreq(data.size, d=timestep)
    
    # And the power (sig_fft is of complex dtype)
    power = np.abs(data)**2
    pos_mask = np.where(sample_freq > 0)
    freqs = sample_freq[pos_mask]
    peak_freq = freqs[power[pos_mask].argmax()]
    
    high_freq_fft = fourier_transform.copy()
    high_freq_fft[np.abs(sample_freq) > peak_freq] = 0
    filtered_sig = fftpack.ifft(high_freq_fft)
    return filtered_sig


def meridional_circulation_index(data):
    """
    Calculates the Meridional Circulation Index (MCI) proposed Francis and Vavrus 2015
    When MCI = 0, the wind is purely zonal, and when MCI= 1 (−1), the flow is from the South (North).
    
           v * abs(v)
    MCI =  ――――――――――
           u**2 * v**2
           
    NOTE: The paper is not clear about whether the absolute value for MCI is taken instead thus 0-1   
    """
    assert 'ua' and 'va' in data.variables, "Cannot compute metric. 'ua' and/or 'va' not found in data" 
    return  data['va']*abs(data['va'])/(data['ua']**2 + data['va']**2)


def get_all_coords_of_jet_occurence(jet_occurence_data):
    all_coords = []
    for val in jet_occurence_data.notnull():
        if val.any():
            for sub_val in val:
                if sub_val:
                    all_coords.append([float(sub_val['lat']), float(sub_val['lon'])])
    return np.array(all_coords)


def get_all_lats_of_jet_centre_for_search(all_coords, lat_resolution, lon_resolution):
    """
        Will get all latitudes that 'may' be a jet-stream centre-point i.e. where latitude appears at least three times for 3*3 lat/lon grid.
        
        NOTE: speeds up calculation as less values are searched through
    """
    ## Step 1. get a counter of all the latitude coordinates
    count_lats = collections.Counter(all_coords[:,0])
    
    ## Step 2. look for all latitudes with at least 3 occurences i.e. X component 3*3 grid
    lats_with_3 = []
    for lat in count_lats.items():
        if lat[1] >= 3:
            lats_with_3.append(lat[0])
    
    ## Step 3. Check if the latitudes above and below the lats with 3 values are present i.e. Y component of for 3*3
    lats_for_search = []
    for lat in lats_with_3:
        if lat-lat_resolution in lats_with_3 and lat+lat_resolution in lats_with_3:
            lats_for_search.append(lat)
    return lats_for_search
            
    
def calculate_jet_centre_points(all_coords, lat_resolution, lon_resolution):
    """
        Will return a list of the coordinates for all jet-centre points
    """

    lats_for_search = get_all_lats_of_jet_centre_for_search(all_coords, lat_resolution, lon_resolution)
    
    jet_centres = []
    for lat in lats_for_search:
        coords_to_search = all_coords[np.where(all_coords[::,0] == lat)]
        for coord in coords_to_search:
            check_for_invalid_coord = coord[0] == 0 or coord[1] == 0 and 360 - lon_resolution not in coords_to_search[::, 1]
            if check_for_invalid_coord:
                continue
            ## check if coord is jet centre point i.e. 9*9 all above 30
            lat_grid_vals = np.arange(coord[0]-lat_resolution, coord[0]+lat_resolution+0.1, lat_resolution)
            lon_grid_vals = np.arange(coord[1]-lon_resolution, coord[1]+lon_resolution+0.1, lon_resolution)
            matrix_vals_to_check = np.array(np.meshgrid(lat_grid_vals, lon_grid_vals)).T.reshape(-1,2)
            matrix_vals_to_check = matrix_vals_to_check % 360 # loop around
            add_coord = True
            for val in matrix_vals_to_check:
                if not val.tolist() in all_coords.tolist():
                    add_coord = False
                    break
            if add_coord:
                jet_centres.append(coord)
            
    return jet_centres
    
    
def get_jet_centre_data(jet_occurence_data):
    """
        Calculates jet-stream centres based on if one jet-stream occurence grid is surrounded by 8 cells of jet-stream occurence (default is 30 m/s)
    """
    ## latitude and longitude resolution
    lat_resolution = float(jet_occurence_data['lat'][1] - jet_occurence_data['lat'][0])
    lon_resolution = float(jet_occurence_data['lon'][1] - jet_occurence_data['lon'][0])
    
    ## make_duplicate data
    jet_centre_data = jet_occurence_data.copy()
    jet_centre_data.loc[:] = np.nan
        
    ## loop over every day
    for time_ind in range(jet_centre_data['time'].size):
        print('On time step: %s' % (time_ind))
        mutliple_time_index = jet_centre_data['time'].size > 1
        
        ## get the lat lon coords of where a jet occurence point is identified
        if mutliple_time_index:
            all_coords = get_all_coords_of_jet_occurence(jet_occurence_data.isel(time=time_ind))
        else:
            all_coords = get_all_coords_of_jet_occurence(jet_occurence_data)
            
        ## calculate jet cores i.e. the central point of a 3*3 grid of jet occurence points
        jet_centres = calculate_jet_centre_points(all_coords, lat_resolution, lon_resolution)

        ## TODO: there's got to be a quicker way
        if mutliple_time_index:
            for centre in jet_centres:
                jet_centre_data.isel(time=time_ind).loc[dict(lat=centre[0], lon=centre[1])] = 1
        else:
            for centre in jet_centres:
                jet_centre_data.loc[dict(lat=centre[0], lon=centre[1])] = 1
        
        ## TODO temporary break if too much data
        if time_ind >= 10:
            print("Stopping calculation as it will take too long, remove this code if more calculation needed!")
            break

    return jet_centre_data