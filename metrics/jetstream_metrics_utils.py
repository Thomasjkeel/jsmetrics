# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to identify or classify jet-stream in the literature.

    This file is in order of publish year of the metrics (see jetstream_metrics_dict.py)
"""

### imports
import numpy as np
import xarray as xr
from scipy import fftpack
import collections
import itertools
from .windspeed_utils import PressureLevelWindSpeedSlice, LatitudeWindSpeedSlice

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def remove_duplicates(vals):
    """
        removes duplicates see: https://stackoverflow.com/questions/2213923/removing-duplicates-from-a-list-of-lists
    """

    vals.sort()
    vals = list(v for v,_ in itertools.groupby(vals))
    return vals


def get_all_plev_hPa(data):
    """
        Will get a list of all the pressure levels in the data in hPa 
    """
    if not 'plev' in data.coords:
        raise KeyError("Data does not contain coord: 'plev'")

    plevs = np.array([plev for plev in data['plev']])
    if data['plev'].units == 'Pa':
        plevs = plevs/100 
    return plevs


def get_sum_weighted_ws(data, all_plevs_hPa):
    """
        Used in Koch et al. 2006
    """
    if not 'plev' in data.coords:
        raise KeyError("Data does not contain coord: 'plev'")

    if not isinstance(all_plevs_hPa, (list, np.ndarray)):
        raise TypeError("array of pressure level needs to be list or numpy.array")

    sum_weighted_ws = 0
    for plev, (i,plev_hPa) in zip(data['plev'], enumerate(all_plevs_hPa)):
        if i != 0:
            plev_hPa = plev_hPa - all_plevs_hPa[i-1]
        sum_weighted_ws += ((data.sel(plev=plev)['ua']** 2 + data.sel(plev=plev)['va']**2)**(1/2)) * plev_hPa
    return sum_weighted_ws


def get_weighted_average_ws(sum_weighted_ws, all_plevs_hPa):
    """
        Used in Koch et al. 2006
    """
    if not isinstance(all_plevs_hPa, (list, np.ndarray)):
        raise TypeError("array of pressure level needs to be list or numpy.array")

    return sum_weighted_ws * (1/(all_plevs_hPa.max() - all_plevs_hPa.min()))


def get_zonal_mean(data):
    """
        Will get the zonal mean either by pressure level (plev) or for one layer
        Used in Woolings et al. 2010
    """
    if not 'lon' in data.coords:
        raise KeyError("data does not contain 'lon' coord")
        
    coords_for_mean = ['lon', 'plev']
    if 'plev' not in data.coords:
        coords_for_mean = ['lon']
    mean_data = data.mean(coords_for_mean)
    return mean_data


def get_latitude_and_speed_where_max_ws(data_row, latitude_col='lat'):
    """
        Will return the latitude and windspeed at the index of maximum wind speed 
        Used in Woolings et al. 2010
    """
    if not data_row.isnull().all():
        max_speed_loc = np.argmax(data_row.data)
        max_speed = data_row[max_speed_loc]
        lat_at_max = float(max_speed[latitude_col].values)
        speed_at_max = float(max_speed.data)
        return lat_at_max, speed_at_max 
    else:
        return None, None


def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.
    
    A low-pass filter removes short-term random fluctations in a time series

    Used in Woolings et al. 2010

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


def apply_lancoz_filter(data, filter_freq, window_size):
    """
        Will carry out Lanczos low-pass filter

        Used in Woolings et al. 2010
    """
    lanczos_weights = low_pass_weights(window_size, 1/filter_freq)
    lanczos_weights_arr = xr.DataArray(lanczos_weights, dims=['window'])
    window_cons = data['ua'].rolling(time=len(lanczos_weights_arr), center=True).construct('window').dot(lanczos_weights_arr)
    return window_cons
    

def fourier_filter(data, timestep=1):
    """
        Carries out a Fourier transform for high frequency filtering
        TAKEN FROM: https://scipy-lectures.org/intro/scipy/auto_examples/plot_fftpack.html
        NOTE: NOT CURRENTLY WORKING PROPERLY
        
        Used in Woolings et al. 2010
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
    peak_freq = freqs[power[pos_mask].argmax()] # TODO change to np.argmax(data) after the tests for Woolings is set up
    
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
        
        For Kuang et al. 2015
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



class JetStreamCoreIdentificationAlgorithm:
    """        
        "As far as object-oriented design is concerned, the breakdown of an algorithm into steps depends on whether the algorithm itself is better seen,
            from its user's point of view, as being executed in one single step, or as multiple steps. If the algorithm is better seen as a single step, 
            the object-oriented implementation can choose to hide all of the implementation details, leaving only the "inputs, process, outputs" visible to the user."
        
        For Manney et al. 2011 etc.
    """
    def __init__(self, data, ws_core_threshold=40, ws_boundary_threshold=30):
        """
            input will need to be longitudinal slice of windspeed values
            
            ws_slice -> one longitude for one day as slice of windspeed
        """
        ## Step 1. make windspeed slice
        self.data = LatitudeWindSpeedSlice(data)
        
        ## Step 2. Get core and potential boundary points
        self.output = self.data.label_slice(self.data['ws'] < ws_core_threshold, 'Core')
        self.output = self.output.where((self.data['ws'] < ws_boundary_threshold) | (self.data['ws'] > ws_core_threshold), other='Potential Boundary')
        
        ## Step 3. Get indexes of jet-stream cores and potential boundaries
        self.core_ids, self.pot_boundary_ids = self.get_indexes_of_core_and_boundaries()
        
        ## Step 4. Set variables needed to keep track of 
        self.current_core_lat = -1
        self.currently_a_core = None

    def __add__(self, other):
        print("TODO: need to implement behaviour for adding slices together")
        pass
        
    
    def __repr__(self):
        """
            Representation of the class. Have it return the labelled data
        """
        print("A total of %d Jet-stream cores have been found in the wind-speed slice" % (len(np.where(self.output['ws'] == 'Core')[1])))
        print("A total of %d potential Jet-stream boundaries have been found in the wind-speed slice" % (len(np.where(self.output['ws'] == 'Potential Boundary')[1])))

        return repr(self.output)
    
    @classmethod
    def run_algorithm(cls, data, ws_core_threshold=40, ws_boundary_threshold=30):
        """
            
        """
        js_algorithm = cls(data, ws_core_threshold=ws_core_threshold, ws_boundary_threshold=ws_boundary_threshold)

        js_core_indexes = js_algorithm.run()
        return js_core_indexes
    
    
    def run(self):
        return self.get_jet_core_boundary()
    
    
    def get_indexes_of_core_and_boundaries(self):
        """
            Will return the indexes in the ws data that ARE jet-stream cores and COULD BE jet-stream core boundaries
        """
        pot_boundary_ids = np.where(self.output['ws'] == 'Potential Boundary')
        core_ids = np.where(self.output['ws'] == 'Core')
        pot_boundary_ids = np.stack(pot_boundary_ids, axis=-1)
        core_ids = np.stack(core_ids, axis=-1)
        return core_ids, pot_boundary_ids

    @staticmethod
    def get_indexes_to_check(pot_boundary):
        """
            Will return an array of indexes to check for potential boundaries or jetstream cores
            Used in Manney et al. 2011
        """
        vals_to_check = []
        if pot_boundary[0] != 0:
            vals_to_check.append([pot_boundary[0]-1, pot_boundary[1]])
        vals_to_check.append([pot_boundary[0]+1, pot_boundary[1]])
        if pot_boundary[1] != 0:
            vals_to_check.append([pot_boundary[0], pot_boundary[1]-1])
        vals_to_check.append([pot_boundary[0], pot_boundary[1]+1])
        return vals_to_check


    def make_pot_jetcore_area(self, vals, area, core_found=False):
        """
        Recursive function that will return the IDs of a jet core boundary i.e. above 30 m/s surrounding a core of 40 m/s
        Will check one an area of potential boundaries contains a core and thus can be called boundaries.

        Used for Manney et al. 2011
        """
        vals_copy = vals.copy()
        for val in vals:
            if val in area:
                continue
            if val in self.core_ids.tolist():
                core_found = True
                ## look for a new core if it is more than 15 degrees away
                if val[1] - self.current_core_lat > 15 and val[1] > self.current_core_lat:
#                     print('THIS IS A NEW CORE')
                    self.current_core_lat = val[1]
                area.append(val)
                new_vals = self.get_indexes_to_check(val)
                vals_copy.extend(new_vals)
                vals_copy = remove_duplicates(vals_copy)
                return self.make_pot_jetcore_area(vals_copy, area=area, core_found=core_found)

            elif val in self.pot_boundary_ids.tolist():
                area.append(val)
                new_vals = self.get_indexes_to_check(val)
                vals_copy.extend(new_vals)
                vals_copy = remove_duplicates(vals_copy)
                return self.make_pot_jetcore_area(vals_copy, area=area, core_found=core_found)
            else:
                vals_copy.remove(val)
                continue
        
        ## reset current core variables
        self.current_core_lat = -1
        self.currently_a_core = None
        return area, core_found


    def get_jet_core_boundary(self):
        """
            Recursive function that will return the IDs of all jet core boundaries i.e. above 30 m/s surrounding a core of 40 m/s
            Will check if an area of potential boundaries contains a core and thus can be called boundaries.

            Used for Manney et al. 2011
        """
        already_covered = []
        js_core_indexes = []
        id_number = 0
        for pot_boundary in self.pot_boundary_ids:
            if pot_boundary.tolist() in already_covered:
                continue
            vals_to_check = self.get_indexes_to_check(pot_boundary)
            area, core_found = self.make_pot_jetcore_area(vals_to_check, area=[])
            already_covered.extend(area)
            already_covered = remove_duplicates(already_covered)
            ## attach area if part of core
            if core_found:
                js_core_indexes.extend([{"id":id_number, "index_of_area":area, "num_of_cores":core_found}])
                id_number += 1


        return js_core_indexes

    
    def work_out_if_two_seperate_cores(self):
        """
            After Manney et al. 2011, will work out if two or more jet cores found within one boundary are seperate or not
        """
        return 


def get_centroid_jet_lat(data, latitude_col='lat'):
    """
        Used in Ceppi et al. 2018
    """
    xs = []
    ys = []
    for lat in data[latitude_col]:
        xs.append(float(lat))
        ys.append(float(data.sel(lat=lat)['ua'].mean()/data['ua'].mean()))
    return np.dot(xs, ys) / np.sum(ys)

