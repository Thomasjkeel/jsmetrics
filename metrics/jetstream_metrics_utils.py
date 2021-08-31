# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to identify or classify jet-stream in the literature.

    This file is in order of publish year of the metrics (see jetstream_metrics_dict.py)
"""

### imports
import numpy as np
import xarray as xr
import scipy.fftpack
import scipy.interpolate
import collections
from .windspeed_utils import PressureLevelWindSpeedSlice, LatitudeWindSpeedSlice
from .general_utils import remove_duplicates, get_num_of_decimal_places, get_local_maxima

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


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
        raise TypeError("array of pressure level needs to be a list or numpy.array")

    return sum_weighted_ws * (1/(all_plevs_hPa.max() - all_plevs_hPa.min()))


def calc_atmospheric_mass_at_kPa(pressure, gravity=9.81, atmospheric_area=5.1e8):
    """
        Will calculate the atmospheric mass at a given pressure level.
        
        Radius of earth (R) = 6.372E3 km;
        Surface area of earth = 4 Pi R^2 = 5.1E8 km^2

        Used in Archer & Caldeira 2008
        
        Returns 

        Parameters
        ---------------
        pressure (float):
            in kPa
        gravity (float):
            m/s^2
    """
    return (pressure/gravity) * atmospheric_area


def get_atm_mass_at_one_hPa(hPa):
    """
        Used in Archer & Caldeira 2008
    """
    kPa = hPa / 10
    atm_mass = calc_atmospheric_mass_at_kPa(kPa) 
    return atm_mass


def get_weighted_average_at_one_Pa(data, Pa, atm_mass):
    """
        Used in Archer & Caldeira 2008
    """
    return atm_mass * (np.sqrt(data['ua'].sel(plev=Pa)**2 + data['va'].sel(plev=Pa)**2))


def get_mass_weighted_average_ws(data, plev_flux=False):
    """
        Used in Archer & Caldeira 2008
    
        TODO: Refactor so neat
    """
    sum_weighted_ws = None #TODO
    for plev_Pa in data['plev'].data:
        plev_hPa = plev_Pa / 100 ## TODO
        atm_mass = get_atm_mass_at_one_hPa(plev_hPa)
        weighted_average = get_weighted_average_at_one_Pa(data, plev_Pa, atm_mass)
        if sum_weighted_ws is None:
            if plev_flux:
                sum_weighted_ws = weighted_average * plev_hPa
            else:
                sum_weighted_ws = weighted_average
        else:
            if plev_flux:
                sum_weighted_ws += weighted_average * plev_hPa
            else:
                sum_weighted_ws += weighted_average
    return sum_weighted_ws


def get_sum_atm_mass(data):
    """
        Used in Archer & Caldeira 2008
    """
    sum_atm_mass = 0
    for plev_Pa in data['plev'].data:
        plev_hPa = plev_Pa / 100 ## TODO
        atm_mass = get_atm_mass_at_one_hPa(plev_hPa)
        sum_atm_mass += atm_mass
    return sum_atm_mass


def calc_mass_weighted_average(data):
    """
        Used in Archer & Caldeira 2008
        TODO: add equation
        TODO: write func desc
    """
    sum_atm_mass = get_sum_atm_mass(data)
    sum_weighted_ws = get_mass_weighted_average_ws(data)
    weighted_average = sum_weighted_ws / sum_atm_mass
    return weighted_average


def calc_mass_flux_weighted_pressure(data):
    """
        Used in Archer & Caldeira 2008
        TODO: add equation
    """
    sum_weighted_ws = get_mass_weighted_average_ws(data)
    sum_weighted_ws_plev_flux = get_mass_weighted_average_ws(data, plev_flux=True)
    mass_flux_weighted_pressure = sum_weighted_ws_plev_flux / sum_weighted_ws
    return mass_flux_weighted_pressure


def calc_mass_flux_weighted_latitude(data, lat_min, lat_max):
    """
        Used in Archer & Caldeira 2008
        TODO: add equation
    """
    assert 'lat' in data.coords, "\'lat\' needs to be in data.coords"
    
    sub_data = data.sel(lat=slice(lat_min, lat_max))
    
    sum_weighted_lat_flux = None
    sum_weighted_ws_by_lat = None
    for lat in sub_data['lat'].data:
        lat_data = sub_data.sel(lat=lat)
        lat_sum_weighted_ws = get_mass_weighted_average_ws(lat_data) 
        if sum_weighted_lat_flux is None:
            sum_weighted_ws_by_lat = lat_sum_weighted_ws
            sum_weighted_lat_flux = lat_sum_weighted_ws * lat
        else:
            sum_weighted_ws_by_lat += lat_sum_weighted_ws
            sum_weighted_lat_flux += lat_sum_weighted_ws * lat
    mass_flux_weighted_latitude = sum_weighted_lat_flux / sum_weighted_ws_by_lat  
    return mass_flux_weighted_latitude


def get_local_jet_maximas_by_day_by_plev(row):
    """
        Used in Schiemann et al 2009
        TODO: add checks
        TODO: will only work is 1 day is the resolution
        TODO: maybe combine with pena-ortiz method
    """
    row['jet_maxima'] = (('plev', 'lat', 'lon'), np.zeros((3, 73, 192))) #TODO
    for lon in row['lon']:
        for plev in row['plev']:
            current = row.sel(lon=lon, plev=plev)
            current = current.where((abs(current['ws']) >= 30) & (current['ua'] > 0))
            local_maxima_lat_inds = get_local_maxima(current['ws'].data)[0]
            if len(local_maxima_lat_inds) > 0:
                for lat_ind in local_maxima_lat_inds:
                    row['jet_maxima'].loc[dict(lat=current['lat'].data[lat_ind], lon=lon, plev=plev)] = 1.0
    return row


def get_zonal_mean(data):
    """
        Will get the zonal mean either by pressure level (plev) or for one layer
        Used in Woolings et al. 2010 & Grise & Polvani 2017
        TODO: add to Archer & Caldiera
    """
    if not 'lon' in data.coords:
        raise KeyError("data does not contain 'lon' coord")
        
    coords_for_mean = ['lon', 'plev']
    if 'plev' not in data.coords or int(data['plev'].count()) == 1:
        coords_for_mean = ['lon']
    zonal_mean = data.mean(coords_for_mean)
    return zonal_mean


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


def apply_lanczos_filter(data, filter_freq, window_size):
    """
        Will carry out Lanczos low-pass filter

        Used in Woolings et al. 2010
    """
    assert filter_freq <= data['time'].count() and filter_freq > 0, "Filter frequency needs to be less\
                                                                     than the number of days in the data\
                                                                    and more than 0 "
    assert window_size <= data['time'].count() and window_size > 0, "Window size needs to be less\
                                                                     than the number of days in the data\
                                                                     and more than 0 "
    assert filter_freq <= window_size, "Filter freq cannot be bigger than window size"

    lanczos_weights = low_pass_weights(window_size, 1/filter_freq)
    lanczos_weights_arr = xr.DataArray(lanczos_weights, dims=['window'])
    window_cons = data['ua'].rolling(time=len(lanczos_weights_arr), center=True).construct('window').dot(lanczos_weights_arr)
    return window_cons
    

def get_latitude_and_speed_where_max_ws(data_row):
    """
        Will return the latitude and windspeed at the index of maximum wind speed from a row of data
        Used in Woolings et al. 2010 & Grise & Polvani 2017
    """
    try:
        assert hasattr(data_row, 'isnull')
    except:
        raise AttributeError("input needs to have isnull method")
    
    if not data_row.isnull().all():
        data_row = data_row.fillna(0.0)
        max_speed_loc = np.argmax(data_row.data)
        max_speed = data_row.isel(lat=max_speed_loc)
        lat_at_max = float(max_speed['lat'].values)
        speed_at_max = float(max_speed.data)
        return lat_at_max, speed_at_max 
    else:
        return None, None



def assign_lat_ws_to_data(data, max_lat_ws):
    """
        Will return a data array with the maximum windspeed and latitude of that 
        maximum wind speed
        Used in Woolings et al. 2010
    """
    max_lats = max_lat_ws[:,0]
    max_ws = max_lat_ws[:,1]
    data_with_max_lats_ws = data.assign({'max_lats':(('time'),max_lats), 'max_ws':(('time'),max_ws)})
    data_with_max_lats_ws['max_lats'] = data_with_max_lats_ws['max_lats'].astype(float)
    data_with_max_lats_ws['max_ws'] = data_with_max_lats_ws['max_ws'].astype(float)
    return data_with_max_lats_ws


def apply_low_freq_fourier_filter(data, highest_freq_to_keep):
    """
        Carries out a Fourier transform for filtering keeping only low frequencies
        ADAPTED FROM: https://scipy-lectures.org/intro/scipy/auto_examples/plot_fftpack.html
        
        Used in Woolings et al. 2010
        Parameters
        ----------
        data : (np.array - 1-d) 
            time series data at regular intervals
        highest_freq_to_keep : (int)
            highest frequency to keep in the fourier transform expression
            NOTE: starts at 0, so highest_freq_to_keep=1 will only keep the constant and first expresion
            
        
        Usage
        ----------
        # Apply filter of the two lowest frequencies
        apply_low_freq_fourier_filter(data, highest_freq_to_keep=2)
            
    """
    ## Fast Fourier Transform on the time series data
    fourier_transform = scipy.fftpack.fft(data)
    
    ## Remove low frequencies
    fourier_transform[highest_freq_to_keep+1:] = 0
    
    ## Inverse Fast Fourier Transform the time series data back
    filtered_sig = scipy.fftpack.ifft(fourier_transform)
    return filtered_sig


def assign_filtered_vals_to_data(data, filtered_max_lats, filtered_max_ws, dim):
    """
        Assigns the filtered data back to the returned dataset
        Used in Woolings et al. 2010
        
    """
    filtered_data = data.assign({'ff_max_lats':((dim), filtered_max_lats),\
                      'ff_max_ws':((dim), filtered_max_ws)})
    filtered_data['ff_max_lats'] = filtered_data['ff_max_lats'].astype(float)
    filtered_data['ff_max_ws'] = filtered_data['ff_max_ws'].astype(float)
    return filtered_data


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
        
        assert ws_core_threshold > ws_boundary_threshold, "core threshold needs to be higher than boundary threshold"
        assert ws_core_threshold > 0, "core threshold needs to be more than 0"
        assert ws_boundary_threshold > 0, "boundary threshold needs to be more than 0"
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


def make_empty_local_wind_maxima_data_var(data):
    """
        Will add a new data var of zeros for local wind maxima
        
        Used in Pena-Ortiz et al. 2013
        
        TODO: add asserts
    """
    data['local_wind_maxima'] = (('time','plev', 'lat', 'lon'), np.zeros((len(data['time']),len(data['plev']),\
                                                                           len(data['lat']),len(data['lon']))))
    return data


def get_potential_local_wind_maximas_by_ws_threshold(ws_slice, ws_threshold):
    """
        Will return a 2-d array of potential local windspeed maximas
        
        Used in Pena-Ortiz et al. 2013
        TODO: add checks
    """
    return ws_slice.where(lambda x: x > 30).fillna(0.0)

    
def get_local_wind_maxima_by_day(row):
    """
        Write function description
        Used in Pena-Ortiz et al. 2013
    """
    
    try:
        assert 'local_wind_maxima' in row.data_vars    
    except Exception as e:
        return print('local_wind_maxima needs to be defined.', e)
    
    for lon in row['lon']:
        current = row.sel(lon=lon)
        pot_local_maximas = get_potential_local_wind_maximas_by_ws_threshold(current['ws'], 30).data
        ind_local_wind_maximas = get_local_maxima(pot_local_maximas, axis=1)
        # Turn into 2-d numpy array 
        ind_local_wind_maximas = np.array([[arr1, arr2] for arr1, arr2 in zip(ind_local_wind_maximas[0], ind_local_wind_maximas[1])])
        for plev_ind, lat_ind in ind_local_wind_maximas:
            row['local_wind_maxima'].loc[dict(lat=current['lat'].data[lat_ind], lon=lon, plev=current['plev'].data[plev_ind])] = 1.0
    return row
          
          
def get_number_of_days_per_monthyear_with_local_wind_maxima(data):
    """
    Will resample by each month and return number of days with local wind maxima
    
    Used in Pena-Ortiz et al. 2013
    """
    data = data['local_wind_maxima'].resample(time='MS').sum().rename({'time':'monthyear'})
    return data


class JetStreamOccurenceAndCentreAlgorithm:
    """
        Have this class inherit from some sort of WS slice of one plev and all Lats + Lons
        
        For Kuang et al. 2015
    """
    
    def __init__(self, data, occurence_ws_threshold=30):
        ## Load in data as a pressure level 2d wind-speed slice
        self.data = PressureLevelWindSpeedSlice(data).values
        assert occurence_ws_threshold > 0, "occurence threshold needs to be more than 0"
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

    
def get_3_latitudes_and_speed_around_max_ws(row):
    """
        Will get the latitude before, on and after where the max windspeed is found
        TODO: think of better name
        TODO: finish func descr
        
        Used in Grise & Polvani 2017 
        
        Parameters
        --------------
        row (xr.DataArray):
    """
    assert 'lat' in row.coords, "\'lat\' needs to be in data.coords"
    
    lat_resolution = float(row['lat'][1] - row['lat'][0])
    lat_min, lat_max = float(row['lat'].min()), float(row['lat'].max())
    max_lat, _ = get_latitude_and_speed_where_max_ws(row)
    neighbouring_lats = get_3_neighbouring_coord_values(max_lat, lat_resolution)
    neighbouring_lats = neighbouring_lats[(neighbouring_lats >= lat_min) & (neighbouring_lats <= lat_max)]
    neighbouring_speeds = row.sel(lat=neighbouring_lats).data 
    return (neighbouring_lats, neighbouring_speeds)


def get_3_neighbouring_coord_values(coord_val, coord_resolution):
    """
        TODO: add to JetStreamOccurenceAndCentreAlgorithm and ...
        
        Used in Grise & Polvani 2017 and ...
        
        Parameters
        --------------
        coord_val (float, int):
            
        coord_resolution (float, int):
            
        Usage
        --------------
        get_3_neighbouring_coord_values(45.0, 1.25)
        >>> [43.75, 45.0, 46.25]
    """
    if type(coord_val) != float or type(coord_resolution) != float:
        coord_val = float(coord_val)
        coord_resolution = float(coord_resolution)
        
    return np.array([coord_val-coord_resolution, coord_val, coord_val+coord_resolution])


def quadratic_func(x, y):
    """
        Used in Grise & Polvani 2017
    """
    p = np.polyfit(x, y, deg=2)
    return p


def apply_quadratic_func(x, y, vals):
    """
        Used in Grise & Polvani 2017
    """
    a, b, c = quadratic_func(x, y)
    return (a * vals**2) + (b * vals) + c


def refine_lat_vals_with_quadratic_func(lats, speeds, lat_vals):
    """
        Will downscale or upscale the resolution of latitude using a quadratic func
        TODO: rename better pls
        
        Used by Grise & Polvani 2017 
    """
    refined_lat_vals = apply_quadratic_func(lats, speeds, lat_vals)
    return refined_lat_vals


def reduce_lat_resolution(lat, resolution):
    """
        Used by Grise & Polvani 2017 & Bracegirdle et al. 2019
    """
    return np.arange(min(lat), max(lat)+resolution, resolution)


def get_latitude_where_max_ws_at_reduced_resolution(lats_and_ws, resolution):
    """
        Makes use of the quadratic func to refine latitude values
        
        Used by Grise & Polvani 2017 
    """ 
    lats, ws = lats_and_ws
    lat_vals =  reduce_lat_resolution(lats, resolution)
    refined_lat_vals = refine_lat_vals_with_quadratic_func(lats, ws, lat_vals)
    decimal_places = get_num_of_decimal_places(resolution)
    return round(lat_vals[np.argmax(refined_lat_vals)], decimal_places)


def get_centroid_jet_lat(data):
    """
        Will get the centroid latitude of the U-wind by day
        
        Only works on data that is 2-d (lat-lon) 
        Used in Ceppi et al. 2018
    """
    xs = []
    ys = []
    for lat in data['lat']:
        xs.append(float(lat))
        ys.append(float(data.sel(lat=lat)['ua'].mean()/data['ua'].mean()))
    return np.dot(xs, ys) / np.sum(ys)


def cubic_spline_interpolation(x, y):
    """
        Used in  Bracegirdle et al. 2019
    """
    return scipy.interpolate.interp1d(x, y, kind='cubic', fill_value='extrapolate')


def run_cubic_spline_interpolation_to_get_max_lat_and_ws(data, resolution, ws_col='ua'):
    """
        Used in  Bracegirdle et al. 2019
        
        Parameters 
        --------------
        data (xr.Dataset):
            must contain coords lat
    """
    refined_lats = reduce_lat_resolution(data['lat'], resolution)
    csi = cubic_spline_interpolation(data['lat'], data[ws_col])
    interpolated_ws = csi(refined_lats)
    max_lat = refined_lats[np.argmax(interpolated_ws)]
    max_ws = max(interpolated_ws)
    return max_lat, max_ws


def run_cubic_spline_interpolation_for_each_climatology_to_get_max_lat_and_ws(data, resolution, time_col):
    max_lats = []
    max_ws = []
    for period in data[time_col].data:
        period_data = data.sel({time_col:period})
        lat, ws = run_cubic_spline_interpolation_to_get_max_lat_and_ws(period_data, resolution=resolution)
        max_lats.append(lat)
        max_ws.append(ws)
    return max_lats, max_ws 