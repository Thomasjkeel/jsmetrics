# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to identify or classify jet-stream in the literature
"""

### imports
import numpy as np
from scipy import fftpack
import collections
import itertools
from .windspeed_utils import PressureLevelWindSpeedSlice, LatitudeWindSpeedSlice

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


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


class JetStreamOccurenceAndCentreAlgorithm:
    """
        Have this class inherit from some sort of WS slice of one plev and all Lats + Lons
    """
    
    def __init__(self, data, occurence_ws_threshold=30):
        ## Load in data as a pressure level 2d wind-speed slice
        self.data = PressureLevelWindSpeedSlice(data).values
        self.jet_occurence = self.data.where(self.data['ws'] >= occurence_ws_threshold)
        self.lat_resolution = float(self.data['lat'][1] - self.data['lat'][0])
        self.lon_resolution = float(self.data['lon'][1] - self.data['lon'][0])
        
        ## make_duplicate data for output
        self.output = self.jet_occurence.copy(deep=True)
        self.output['ws'].loc[:] = np.nan
        
        ## Initialise lists needed for search algorithm
        self.all_coords = []
        self.lats_with_3 = []
        self.lats_for_search = []
        self.jet_centres = []
        
    @classmethod
    def run_algorithm(cls, data):
        return cls(data).run()
    
    def run(self):
        print('starting')
        self.get_all_coords_of_jet_occurence()
        self.all_coords_arr = np.array(self.all_coords)
        ## Get a counter of all the latitude coordinates
        self.count_lats = collections.Counter(self.all_coords_arr[:,0])
        self.get_all_lats_of_jet_centre_for_search()
        self.calculate_jet_centre_points()
        self.get_jet_centre_data()
        print('done!')
    
    def get_jet_centre_data(self):
        """
            Calculates jet-stream centres based on if one jet-stream occurence grid is surrounded by 8 cells of jet-stream occurence (default is 30 m/s)
        """
        ## TODO: there's got to be a quicker way
        for centre in self.jet_centres:
            self.output['ws'].loc[dict(lat=centre[0], lon=centre[1])] = 1
    
    
    def get_all_coords_of_jet_occurence(self):
        for val in self.jet_occurence['ws'].notnull():
            if val.any():
                for sub_val in val:
                    if sub_val:
                        self.all_coords.append([float(sub_val['lat']), float(sub_val['lon'])])

    def get_all_lats_of_jet_centre_for_search(self):
        """
            Will get all latitudes that 'may' be a jet-stream centre-point i.e. where latitude appears at least three times for 3*3 lat/lon grid.
            NOTE: speeds up calculation as less values are searched through
        """
        ## Step 1. look for all latitudes with at least 3 occurences i.e. X component 3*3 grid
        self.get_all_latitudes_that_occur_at_least_three_times()
        ## Step 2. Check if the latitudes above and below the lats with 3 values are present i.e. Y component of for 3*3
        self.get_all_latitudes_available_in_3by3_grid()
    
    def get_all_latitudes_that_occur_at_least_three_times(self):
        for lat in self.count_lats.items():
            if lat[1] >= 3:
                self.lats_with_3.append(lat[0])
    
    def get_all_latitudes_available_in_3by3_grid(self):
        for lat in self.lats_with_3:
            if lat-self.lat_resolution in self.lats_with_3 and lat+self.lat_resolution in self.lats_with_3:
                self.lats_for_search.append(lat)

    def calculate_jet_centre_points(self):
        """
            Will return a list of the coordinates for all jet-centre points
        """
        for lat in self.lats_for_search:
            coords_to_search = self.all_coords_arr[np.where(self.all_coords_arr[::,0] == lat)]
            for coord in coords_to_search:
                if coord[0] == 0 or coord[1] == 0 and 360 - self.lon_resolution not in coords_to_search[::, 1]:
                    continue
                ## check if coord is jet centre point i.e. 9*9 all above 30
                lat_grid_vals = np.arange(coord[0]-self.lat_resolution, coord[0]+self.lat_resolution+0.1, self.lat_resolution)
                lon_grid_vals = np.arange(coord[1]-self.lon_resolution, coord[1]+self.lon_resolution+0.1, self.lon_resolution)
                matrix_vals_to_check = np.array(np.meshgrid(lat_grid_vals, lon_grid_vals)).T.reshape(-1,2)
                matrix_vals_to_check = matrix_vals_to_check % 360 # loop around
                add_coord = True
                for val in matrix_vals_to_check:
                    if not val.tolist() in self.all_coords:
                        add_coord = False
                        break
                if add_coord:
                    self.jet_centres.append(coord)
