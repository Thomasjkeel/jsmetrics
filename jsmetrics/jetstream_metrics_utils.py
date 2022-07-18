# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to
    identify or classify jet-stream in the literature.

    This file is in order of publish year of the metrics
    (see details_for_all_metrics.py)
"""

# imports
import collections
import numpy as np
import matplotlib.pyplot
import xarray as xr
import scipy.fftpack
import scipy.interpolate
import shapely.geometry
from . import windspeed_utils, general_utils

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def get_sum_weighted_ws(data, all_plevs_hPa):
    """
    Get sum of weighted windspeed.
    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    sum weighted windspeed = integral(p2, p1)(u^2+v^2)^(1/2)dp
    where p1, p2 is min, max pressure level

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind
    all_plevs_hPa : array-like
        list of hPa unit pressure levels

    Returns
    ----------
    sum_weighted_ws : xarray.Dataset
        Data containing sum weighted windspeed values
    """
    if "plev" not in data.coords:
        raise KeyError("Data does not contain coord: 'plev'")

    if not isinstance(all_plevs_hPa, (list, np.ndarray)):
        raise TypeError(
            "array of pressure level needs to be list or numpy.array"
        )

    sum_weighted_ws = 0
    for plev, (i, plev_hPa) in zip(data["plev"], enumerate(all_plevs_hPa)):
        if i != 0:
            plev_hPa = plev_hPa - all_plevs_hPa[i - 1]
        sum_weighted_ws += (
            (data.sel(plev=plev)["ua"] ** 2 + data.sel(plev=plev)["va"] ** 2)
            ** (1 / 2)
        ) * plev_hPa
    return sum_weighted_ws


def get_weighted_average_ws(sum_weighted_ws, all_plevs_hPa):
    """
    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    weighted average windspeed = 1/(p2-p1) * sum average windspeed
    where p1, p2 is min, max pressure level

    Parameters
    ----------
    sum_weighted_ws : xarray.Dataset
        Data containing sum weighted windspeed values
    all_plevs_hPa : array-like
        list of hPa unit pressure levels

    Returns
    ----------
    weighted_average_ws : xarray.Dataset
        Data containing weighted average windspeed values
    """
    if not isinstance(all_plevs_hPa, (list, np.ndarray)):
        raise TypeError(
            "array of pressure level needs to be a list or numpy.array"
        )
    weighted_average_ws = sum_weighted_ws * (
        1 / (all_plevs_hPa.max() - all_plevs_hPa.min())
    )
    return weighted_average_ws


def calc_atmospheric_mass_at_kPa(
    pressure, gravity=9.81, atmospheric_area=5.1e8
):
    """
    Will calculate the atmospheric mass at a given pressure level.

    Radius of earth (R) = 6.372E3 km;
    Surface area of earth = 4 Pi R^2 = 5.1E8 km^2

    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Returns

    Parameters
    ---------------
    pressure : float
        in kPa
    gravity : float
        m/s^2
    """
    return (pressure / gravity) * atmospheric_area


def get_atm_mass_at_one_hPa(hPa):
    """
    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Parameters
    ----------
    hPa : int or float
        One pressure level in hPa

    Returns
    ----------
    atm_mass : int or float
        Atmospheric mass at given hPa pressure level
    """
    kPa = hPa / 10
    atm_mass = calc_atmospheric_mass_at_kPa(kPa)
    return atm_mass


def get_weighted_average_at_one_Pa(data, Pa, atm_mass, ws_col):
    """
    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Parameters
    ----------
    data : xarray.Dataset
        Data containing windspeed ('ws') or u-, v-component wind
    Pa : int or float
        Pressure level in Pascals
    atm_mass : int or float
        Atmospheric mass at given hPa pressure level
    ws_col : string
        Name of column to calculate weighted average from (e.g. 'ws', 'ua', 'va')

    Returns
    ----------
    output : xarray.Dataset
        Data with weighted average at a single pressure level
    """
    return atm_mass * (data[ws_col].sel(plev=Pa))


def get_mass_weighted_average_wind(data, ws_col, plev_flux=False):
    """
    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Get mass-weighted average wind-speed from 'u', 'v' component wind or 'wind vector'.

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind
    ws_col : string
        Name of column to calculate weighted average from (e.g. 'ws', 'ua', 'va')

    Returns
    ----------
    output : xarray.Dataset
        Data containing mass weighted average ws

    Usage
    ----------
        mon_mean = data.groupby("time.month").mean()
        mass_weighted_average = jetstream_metrics_utils.get_mass_weighted_average_ws(mon_mean)
    """
    general_utils.check_at_least_two_plevs_in_data(data)
    sum_weighted_ws = None  # TODO
    for plev_Pa in data["plev"].data:
        plev_hPa = plev_Pa / 100  # TODO
        atm_mass = get_atm_mass_at_one_hPa(plev_hPa)
        weighted_average = get_weighted_average_at_one_Pa(
            data, plev_Pa, atm_mass, ws_col
        )
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
    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Parameters
    ----------
    data : xarray.Dataset
        Data containing plev

    Returns
    ----------
    sum_atm_mass : int or float
        Sum of atmospheric mass
    """
    sum_atm_mass = 0
    general_utils.check_at_least_two_plevs_in_data(data)
    for plev_Pa in data["plev"].data:
        plev_hPa = plev_Pa / 100  # TODO
        atm_mass = get_atm_mass_at_one_hPa(plev_hPa)
        sum_atm_mass += atm_mass
    return sum_atm_mass


def calc_mass_weighted_average(data, ws_col):
    """
    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    TODO: add equation

    Parameters
    ----------
    data : xarray.Dataset
        Data containing plev
    ws_col : string
        Name of column to calculate weighted average from (e.g. 'ws', 'ua', 'va')

    Returns
    ----------
    weighted_average : xr.DataArray
        Data with weighted average windspeed based on sum atmospheric mass
    """
    sum_atm_mass = get_sum_atm_mass(data)
    sum_weighted_ws = get_mass_weighted_average_wind(data, ws_col)
    weighted_average = sum_weighted_ws / sum_atm_mass
    return weighted_average


def calc_mass_flux_weighted_pressure(data, ws_col):
    """
    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    TODO: add equation

    Parameters
    ----------
    data : xarray.Dataset
        Data containing windspeed and plev
    ws_col : string
        Name of column to calculate weighted average from (e.g. 'ws', 'ua', 'va')

    Returns
    ----------
    mass_flux_weighted_pressure : xr.DataArray
        Data with mass flux weighted pressure
    """
    sum_weighted_ws = get_mass_weighted_average_wind(data, ws_col=ws_col)
    sum_weighted_ws_plev_flux = get_mass_weighted_average_wind(
        data, ws_col=ws_col, plev_flux=True
    )
    mass_flux_weighted_pressure = sum_weighted_ws_plev_flux / sum_weighted_ws
    return mass_flux_weighted_pressure


def calc_mass_flux_weighted_latitude(data, lat_min, lat_max, ws_col):
    """
    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    TODO: add equation
    WARNING: Problem with including 1000 hPa

    Parameters
    ----------
    data : xarray.Dataset
        Data containing windspeed and lat
    lat_min : int or float
        Minimum latitude to consider for weighted latitude
    lat_max : int or float
        Maximum latitude to consider for weighted latitude
    ws_col : string
        Name of column to calculate weighted average from (e.g. 'ws', 'ua', 'va')

    Returns
    ----------
    mass_flux_weighted_latitude : xr.DataArray
        Data with mass flux weighted latitude
    """
    assert "lat" in data.coords, "'lat' needs to be in data.coords"

    sub_data = data.sel(lat=slice(lat_min, lat_max))

    sum_weighted_lat_flux = None
    sum_weighted_ws_by_lat = None
    for lat in sub_data["lat"].data:
        lat_data = sub_data.sel(lat=lat)
        lat_sum_weighted_ws = get_mass_weighted_average_wind(lat_data, ws_col)
        if sum_weighted_lat_flux is None:
            sum_weighted_ws_by_lat = lat_sum_weighted_ws
            sum_weighted_lat_flux = lat_sum_weighted_ws * lat
        else:
            sum_weighted_ws_by_lat += lat_sum_weighted_ws
            sum_weighted_lat_flux += lat_sum_weighted_ws * lat
    mass_flux_weighted_latitude = (
        sum_weighted_lat_flux / sum_weighted_ws_by_lat
    )
    return mass_flux_weighted_latitude


def get_local_jet_maximas_by_timeunit_by_plev(row):
    """
    Component of method from Schiemann et al 2009 https://doi.org/10.1175/2008JCLI2625.1
    TODO: add checks
    TODO: will only work is 1 day is the resolution
    TODO: maybe combine with pena-ortiz method

    Parameters
    ----------
    row : xarray.Dataset
        Data of a sinlge time unit containing windspeed (ws), plev, lat, lon

    Returns
    ----------
    row : xarray.Dataset
        Data of a sinlge time unit with value for jet-maxima (1 == maxima, 0 == none)

    """
    row["jet_maxima"] = (
        ("plev", "lat", "lon"),
        np.zeros((row["plev"].size, row["lat"].size, row["lon"].size)),
    )  # TODO
    for lon in row["lon"]:
        for plev in row["plev"]:
            current = row.sel(lon=lon, plev=plev)
            current = current.where(
                (abs(current["ws"]) >= 30) & (current["ua"] > 0)
            )
            local_maxima_lat_inds = general_utils.get_local_maxima(
                current["ws"].data
            )[0]
            if len(local_maxima_lat_inds) > 0:
                for lat_ind in local_maxima_lat_inds:
                    row["jet_maxima"].loc[
                        dict(
                            lat=current["lat"].data[lat_ind],
                            lon=lon,
                            plev=plev,
                        )
                    ] = 1.0
    return row


def get_zonal_mean(data):
    """
    Component of method from Woolings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Barnes & Polvani (2015) http://journals.ametsoc.org/doi/10.1175/JCLI-D-14-00589.1
    & Grise & Polvani (2017) https://doi.org/10.1175/JCLI-D-16-0849.1
    & Ceppi et al (2018) https://doi.org/10.1175/JCLI-D-17-0323.1

    Will get the zonal mean either by pressure level (plev) or for one layer
    TODO: add to Archer & Caldiera

    Parameters
    ----------
    data : xarray.Dataset
        Data containing lon and plev coords

    Returns
    ----------
    zonal_mean : xarray.DataSet
        zonal mean data

    Raises
    ----------
    KeyError
        when 'lon' not discovered as coord
    """
    if "lon" not in data.coords:
        raise KeyError("data does not contain 'lon' coord")

    coords_for_mean = ["lon", "plev"]
    if "plev" not in data.coords or int(data["plev"].count()) == 1:
        coords_for_mean = ["lon"]
    zonal_mean = data.mean(coords_for_mean)
    return zonal_mean


def calc_low_pass_weights(window, cutoff):
    """
    Component of method from Woolings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Calculate weights for a low pass Lanczos filter.
    A low-pass filter removes short-term random fluctations in a time series

    TAKEN FROM:
    https://scitools.org.uk/iris/docs/v1.2/examples/graphics/SOI_filtering.html

    Parameters
    ----------
    window : int
        The length of the filter window.
    cutoff : float
        The cutoff frequency in inverse time steps.

    Returns
    ----------
    output : numpy.array
        filtered weights
    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1.0, n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2.0 * np.pi * cutoff * k) / (np.pi * k)
    w[n - 1 : 0 : -1] = firstfactor * sigma
    w[n + 1 : -1] = firstfactor * sigma
    return w[0 + (window % 2) : -1]  # edited from w[1:-1]


def apply_lanczos_filter(dataarray, filter_freq, window_size):
    """
    Component of method from Woolings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Will carry out Lanczos low-pass filter

    Parameters
    ----------
    datarray : xarray.DataArray
        Data to apply filter to containing zonal mean (u-component) wind
    filter_freq : int
        number of days in filter
    window_size : int
        number of days in window for Lancoz filter

    Returns
    ----------
    window_cons : xarray.DataArray
        Filtered zonal mean data
    """
    if (
        dataarray["time"].count() <= filter_freq
        or dataarray["time"].count() <= window_size
        or window_size <= filter_freq
    ):
        raise ValueError(
            "Time series is too short to apply %s window for Lanczos filter"
            % (window_size)
        )

    assert (
        filter_freq >= 0 and window_size >= 0
    ), "both filter_freq and window need to be more than 0"
    assert isinstance(
        dataarray, xr.DataArray
    ), "Input data needs to be a data array"

    lanczos_weights = calc_low_pass_weights(window_size, 1 / filter_freq)
    lanczos_weights_arr = xr.DataArray(lanczos_weights, dims=["window"])
    window_cons = (
        dataarray.rolling(time=len(lanczos_weights_arr), center=True)
        .construct("window")
        .dot(lanczos_weights_arr)
    )
    return window_cons


def get_latitude_and_speed_where_max_ws(data_row):
    """
    Component of method from Woolings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    Barnes & Simpson 2017 https://doi.org/10.1175/JCLI-D-17-0299.1
    & Grise & Polvani 2017 https://doi.org/10.1175/JCLI-D-16-0849.1

    Will return the latitude and windspeed at the index of maximum wind speed
    from a row of data
    TODO: expects only one cell to be max

    Parameters
    ----------
    data_row : xarray.DataArray
        Data of single time unit containing u-component wind

    Returns
    ----------
    lat_at_max : int or float
        Latitude of maximum windspeed
    speed_at_max : int or float
        Speed at latitude of maximum windspeed
    """
    try:
        assert hasattr(data_row, "isnull")
    except Exception as e:
        raise AttributeError("input needs to have isnull method.") from e

    if not data_row.isnull().all():
        data_row = data_row.fillna(0.0)
        max_ws = data_row.where(
            np.abs(data_row) == np.abs(data_row).max(), drop=True
        ).squeeze()
        if max_ws.size > 1:
            print("Warning: more than one max value found, picking the max!")
            max_ws = max_ws[0].max()
        lat_at_max = float(max_ws["lat"])
        speed_at_max = float(max_ws.data)
        return lat_at_max, speed_at_max
    else:
        return np.nan, np.nan


def assign_jet_lat_and_speed_to_data(
    data, max_lat_ws, max_lats_col="jet_lat", max_ws_col="jet_speed"
):
    """
    Component of method from Woolings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Simpson 2017 https://doi.org/10.1175/JCLI-D-17-0299.1

    Will return a data array with the maximum windspeed and latitude of that maximum wind speed

    Parameters
    ----------
    data : xarray.DataArray
        Data of single time unit containing u-component wind
    max_lat_ws : np.array
        Array of maximum latitude and windspeed at those maximums for each timeunit

    Returns
    ----------
    data_with_max_lats_ws : xarray.Dataset
        Data with maximum latitude and windspeed attached
    """
    max_lats = max_lat_ws[:, 0]
    max_ws = max_lat_ws[:, 1]
    data_with_max_lats_ws = data.assign(
        {max_lats_col: (("time"), max_lats), max_ws_col: (("time"), max_ws)}
    )
    data_with_max_lats_ws[max_lats_col] = data_with_max_lats_ws[
        max_lats_col
    ].astype(float)
    data_with_max_lats_ws[max_ws_col] = data_with_max_lats_ws[
        max_ws_col
    ].astype(float)
    return data_with_max_lats_ws


def apply_low_freq_fourier_filter(data, highest_freq_to_keep):
    """
    Component of method from Woolings et al (2010) http://dx.doi.org/10.1002/qj.625

    Carries out a Fourier transform for filtering keeping only low frequencies
    ADAPTED FROM:
    https://scipy-lectures.org/intro/scipy/auto_examples/plot_fftpack.html

    Parameters
    ----------
    data : np.array - 1-d
        time series data at regular intervals
    highest_freq_to_keep : int
        highest frequency to keep in the fourier transform expression
        NOTE: starts at 0, so highest_freq_to_keep=1
              will only keep the constant and first expresion

    Returns
    ----------
    filtered_sig : np.array
        Fourier-filtered signal

    Usage
    ----------
    # Apply filter of the two lowest frequencies
    apply_low_freq_fourier_filter(data, highest_freq_to_keep=2)
    """
    # Fast Fourier Transform on the time series data
    fourier_transform = scipy.fftpack.fft(data)

    # Remove low frequencies
    fourier_transform[highest_freq_to_keep + 1 :] = 0

    # Inverse Fast Fourier Transform the time series data back
    filtered_sig = scipy.fftpack.ifft(fourier_transform)
    return filtered_sig


def assign_filtered_lats_and_ws_to_data(
    data, filtered_max_lats, filtered_max_ws, dim
):
    """
    Component of method from Woolings et al (2010) http://dx.doi.org/10.1002/qj.625
    Assigns the filtered data back to the returned dataset

    Parameters
    ----------
    data : xarray.Dataset
        Data to have fourier filtered data assigned to
    filtered_max_lats : numpy.array
        Fourier-filtered maximum latitude by given timeunit
    filtered_max_lats : numpy.array
        Fourier-filtered maximum speed at maximum latitude by given timeunit

    Returns
    ----------
    filtered_data : xarray.Dataset
        Data with fourier-filtered maximum lat and windspeed attached to it
    """
    filtered_data = data.assign(
        {
            "ff_jet_lat": ((dim), filtered_max_lats),
            "ff_jet_speed": ((dim), filtered_max_ws),
        }
    )
    filtered_data["ff_jet_lat"] = filtered_data["ff_jet_lat"].astype(float)
    filtered_data["ff_jet_speed"] = filtered_data["ff_jet_speed"].astype(float)
    return filtered_data


def run_jet_core_algorithm_on_one_day(
    row, ws_core_threshold, ws_boundary_threshold
):
    """
    Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
    Runs JetStreamCoreIdentificationAlgorithm method on a single time unit

    Parameters
    ----------
    row : xarray.Dataset
        Data of single time unit containing u- and v-component wind
    ws_core_threshold : int or float
        Threshold used for jet-stream core point
    ws_boundary_threshold : int or float
        Threshold for jet-stream boundary point

    Returns
    ----------
    row : xarray.Dataset
        Data for one time unit containing jet-cores (ID number relates to each unique core)
    """
    row["jet_core_id"] = (
        ("plev", "lat", "lon"),
        np.zeros((row["plev"].size, row["lat"].size, row["lon"].size)),
    )  # TODO
    for lon in row["lon"]:
        current = row.sel(lon=lon)
        core_alg = JetStreamCoreIdentificationAlgorithm(
            current,
            ws_core_threshold=ws_core_threshold,
            ws_boundary_threshold=ws_boundary_threshold,
        )
        core_alg.run()
        row["jet_core_id"].loc[dict(lon=lon)] = core_alg.output_data["core_id"]
    return row


class JetStreamCoreIdentificationAlgorithm:
    """
    Jet-stream core identification algorithm
    Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
    """

    def __init__(self, data, ws_core_threshold=40, ws_boundary_threshold=30):
        """
        Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
        input will need to be longitudinal slice of windspeed values


        Parameters
        ----------
        data : xarray.Dataset
            Data of single time unit containing u- and v-component wind
        ws_core_threshold : int or float
            Threshold used for jet-stream core point (default=40)
        ws_boundary_threshold : int or float
            Threshold for jet-stream boundary point (default=30)
        """
        try:
            assert (
                ws_core_threshold > ws_boundary_threshold
                and ws_core_threshold > 0
                and ws_boundary_threshold > 0
            )
        except Exception as e:
            raise ValueError(
                "Windspeed core threshold needs to be more than boundary\
                    threshold and both need to be more than 0"
            ) from e
        # standardise data
        data = general_utils.standardise_dimension_order(
            data, dim_order=(..., "lat", "plev")
        )
        # Step 1. make windspeed slice
        self._lat_ws_slice = windspeed_utils.LatitudeWindSpeedSlice(data)
        # Step 2. Get core and potential boundary points
        self._labelled_data = self._lat_ws_slice.label_slice(
            self._lat_ws_slice["ws"] < ws_core_threshold, "Core"
        )
        self._labelled_data = self._labelled_data.where(
            (self._lat_ws_slice["ws"] < ws_boundary_threshold)
            | (self._lat_ws_slice["ws"] > ws_core_threshold),
            other="Potential Boundary",
        )

        # Step 3. Get indexes of jet-stream cores and potential boundaries
        (
            self._initial_core_ids,
            self._pot_boundary_ids,
        ) = self._get_indexes_of_core_and_boundaries()

        # Step 4. Set variables needed to keep track of
        self._current_core_lat = -1
        self._currently_a_core = None
        self.algorithm_has_run = False

    def __add__(self, other):
        print("TODO: need to implement behaviour for adding slices together")
        pass

    def __repr__(self):
        """
        Representation of the class. Have it return the labelled data
        """
        if not self.algorithm_has_run:
            print(
                "A total of %d initial Jet-stream cores have been found\
                 in the wind-speed slice"
                % (
                    self._labelled_data["ws"]
                    .where(lambda x: x == "Core")
                    .count()
                )
            )
            print(
                "A total of %d potential Jet-stream boundaries have been found\
                 in the wind-speed slice"
                % (
                    self._labelled_data["ws"]
                    .where(lambda x: x == "Potential Boundary")
                    .count()
                )
            )
            return repr(self._labelled_data)
        else:
            print(
                "A total of %d cores have been discovered"
                % (self.num_of_cores)
            )
            return repr(self.output_data)

    @classmethod
    def run_algorithm(
        cls, data, ws_core_threshold=40, ws_boundary_threshold=30
    ):
        """
        Class method for running algorithm
        """
        js_algorithm = cls(
            data,
            ws_core_threshold=ws_core_threshold,
            ws_boundary_threshold=ws_boundary_threshold,
        )

        js_core_indexes = js_algorithm.run()
        return js_core_indexes

    def run(self):
        self.final_jet_cores = self._get_jet_core_boundary()
        self.output_data = self._add_cores_to_data()
        self.algorithm_has_run = True

    def _get_indexes_of_core_and_boundaries(self):
        """
        Will return the indexes in the ws data that ARE jet-stream cores
        and COULD BE jet-stream core boundaries
        """
        pot_boundary_ids = np.where(
            self._labelled_data["ws"] == "Potential Boundary"
        )
        initial_core_ids = np.where(self._labelled_data["ws"] == "Core")
        pot_boundary_ids = np.stack(pot_boundary_ids, axis=-1)
        initial_core_ids = np.stack(initial_core_ids, axis=-1)
        return initial_core_ids, pot_boundary_ids

    @staticmethod
    def _get_indexes_to_check(pot_boundary):
        """
        Will return an array of indexes to check for potential boundaries
        or jetstream cores.
        Used in Manney et al. 2011
        """
        vals_to_check = []
        if pot_boundary[0] != 0:
            vals_to_check.append([pot_boundary[0] - 1, pot_boundary[1]])
        vals_to_check.append([pot_boundary[0] + 1, pot_boundary[1]])
        if pot_boundary[1] != 0:
            vals_to_check.append([pot_boundary[0], pot_boundary[1] - 1])
        vals_to_check.append([pot_boundary[0], pot_boundary[1] + 1])
        return vals_to_check

    def _get_pot_jetcore_area(self, vals, area, core_found=False):
        """
        Recursive function that will return the IDs of a jet core boundary
        i.e. above 30 m/s surrounding a core of 40 m/s
        Will check one an area of potential boundaries contains a core
        and thus can be called boundaries.
        """
        vals_copy = vals.copy()
        for val in vals:
            if val in area:
                continue
            if val in self._initial_core_ids.tolist():
                core_found = True
                # look for a new core if it is more than 15 degrees away
                if (
                    val[1] - self._current_core_lat > 15
                    and val[1] > self._current_core_lat
                ):
                    #  print('THIS IS A NEW CORE')
                    self._current_core_lat = val[1]
                area.append(val)
                new_vals = self._get_indexes_to_check(val)
                vals_copy.extend(new_vals)
                vals_copy = general_utils.remove_duplicates(vals_copy)
                return self._get_pot_jetcore_area(
                    vals_copy, area=area, core_found=core_found
                )

            elif val in self._pot_boundary_ids.tolist():
                area.append(val)
                new_vals = self._get_indexes_to_check(val)
                vals_copy.extend(new_vals)
                vals_copy = general_utils.remove_duplicates(vals_copy)
                return self._get_pot_jetcore_area(
                    vals_copy, area=area, core_found=core_found
                )
            else:
                vals_copy.remove(val)
                continue

        # reset current core variables
        self._current_core_lat = -1
        self._currently_a_core = None
        return area, core_found

    def _get_jet_core_boundary(self):
        """
        Recursive function that will return the IDs of all jet core boundaries
        i.e. above 30 m/s surrounding a core of 40 m/s
        Will check if an area of potential boundaries contains a core
        and thus can be called boundaries.
        """
        already_covered = []
        js_core_indexes = []
        id_number = 0
        for pot_boundary in self._pot_boundary_ids:
            if pot_boundary.tolist() in already_covered:
                continue
            vals_to_check = self._get_indexes_to_check(pot_boundary)
            area, core_found = self._get_pot_jetcore_area(
                vals_to_check, area=[]
            )
            already_covered.extend(area)
            already_covered = general_utils.remove_duplicates(already_covered)
            # add area to js_core_indexes if part of core
            if core_found:
                id_number += 1
                js_core_indexes.extend(
                    [{"id": id_number, "index_of_area": area}]
                )
        self.num_of_cores = id_number
        return js_core_indexes

    def _add_cores_to_data(self):
        self._lat_ws_slice.values["core_id"] = (
            ("plev", "lat"),
            np.zeros(
                (
                    self._lat_ws_slice.values["plev"].size,
                    self._lat_ws_slice.values["lat"].size,
                )
            ),
        )
        for jet_core in self.final_jet_cores:
            for lat, plev in jet_core["index_of_area"]:
                self._lat_ws_slice.values["core_id"].loc[
                    dict(
                        lat=self._lat_ws_slice.values["lat"].data[lat],
                        plev=self._lat_ws_slice.values["plev"].data[plev],
                    )
                ] = jet_core["id"]
        return self._lat_ws_slice.values


def calc_jet_width_for_one_day(data_row, jet_lat, jet_speed):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Calculates jet width on a given data array using the calculated jet lat and speed at that lat

    Parameters
    ----------
    data_row : xarray.DataArray
        Data of single time unit containing u-component wind
    jet_lat : int or float
        Latitude of jet
    jet_speed : int or float
        Speed of jet at jet_lat

    Returns
    ----------
    jet_widths : array-like
        Width of jet-stream in latitude. Units for latitude are the same as the input units i.e. likely to be 'degrees_north'
    """
    lat_resolution = float(data_row["lat"][1] - data_row["lat"][0])
    if not jet_speed:
        return np.nan
    possible_surrounding_jet_vals = (
        get_possible_surrounding_jet_lat_vals_for_one_day(data_row, jet_speed)
    )
    if possible_surrounding_jet_vals.size == data_row["lat"].size:
        print("No jet-width determined")
        return np.nan
    possible_surrounding_jet_lat_vals = possible_surrounding_jet_vals["lat"]
    jet_width = get_width_of_jet_from_surrounding_lat_vals(
        jet_lat, possible_surrounding_jet_lat_vals, lat_resolution
    )
    return jet_width


def get_possible_surrounding_jet_lat_vals_for_one_day(
    one_day_wind_data, jet_speed
):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Returns array of windspeed values that are within half speed of the jet_speed given

    Parameters
    ----------
    one_day_wind_data : xarray.DataArray
        Data of single time unit containing u-component wind
    jet_speed : int or float
        Speed of jet at jet_lat

    Returns
    ----------
    possible_surrounding_jet_vals : array-like
        DataArray of True or False as to whether the windspeed value is within half speed of the jet_speed
    """
    half_jet_speed = jet_speed / 2
    possible_surrounding_jet_vals = one_day_wind_data[
        one_day_wind_data > half_jet_speed
    ]
    return possible_surrounding_jet_vals


def get_width_of_jet_from_surrounding_lat_vals(
    jet_lat, possible_surrounding_jet_lat_vals, lat_resolution
):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Get latitude width of jet around the given jet latitude value. Width given in units of latitude

    Parameters
    ----------
    jet_lat : int or float
        Latitude of jet
    possible_surrounding_jet_lat_vals : array-like
        DataArray of lats where the windspeed value is within half speed of the jet_speed
    lat_resolution : int or float
        Latitude resolution in degrees

    Returns
    ----------
    jet_width : array-like
        Width of jet-stream in latitude. Unit for latitude are the same as the input units i.e. likely to be 'degrees_north'
    """
    lats_within_jet = get_all_surrounding_jet_lats_around_max_lat(
        jet_lat, possible_surrounding_jet_lat_vals, lat_resolution
    )
    if not lats_within_jet:
        return np.nan
    jet_width = max(lats_within_jet) - min(lats_within_jet)
    return jet_width


def get_all_surrounding_jet_lats_around_max_lat(
    jet_lat, possible_surrounding_jet_lat_vals, lat_resolution
):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Returns all the latitudes in a continuous group (i.e. within one unit of lat resolution) around jet lat

    Parameters
    ----------
    jet_lat : int or float
        Latitude of jet
    possible_surrounding_jet_lat_vals : array-like
        DataArray of lats where the windspeed value is within half speed of the jet_speed
    lat_resolution : int or float
        Latitude resolution in degrees

    Returns
    ----------
    lats_within_jet : list
        All latitudes in continous group (i.e. within lat_resolution) around jet_lat
    """
    lats_within_jet = []
    lats_at_higher_lats = check_lats_within_jet_in_one_direction(
        jet_lat, possible_surrounding_jet_lat_vals, lat_resolution
    )
    lats_at_lower_lats = check_lats_within_jet_in_one_direction(
        jet_lat, possible_surrounding_jet_lat_vals, -lat_resolution
    )
    lats_within_jet.extend(lats_at_higher_lats)
    lats_within_jet.extend(lats_at_lower_lats)
    return lats_within_jet


def check_lats_within_jet_in_one_direction(
    jet_lat, possible_surrounding_jet_lat_vals, amount_to_add_to_lat_val
):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Creates a list of latitudes that are within a given distance from the jet_lat and only in one direction i.e. increasing or decreasing

    Parameters
    ----------
    jet_lat : int or float
        Latitude of jet
    possible_surrounding_jet_lat_vals : array-like
        DataArray of lats where the windspeed value is within half speed of the jet_speed
    amount_to_add_to_lat_val : int or float
        Check if lat within this amount

    Returns
    ----------
    lats_within_jet : list
        All latitudes in continous group (i.e. within lat_resolution) around jet_lat
    """
    lats_within_jet = []
    current_lat = jet_lat
    current_lat_is_within_range = True
    while current_lat_is_within_range:
        if current_lat in possible_surrounding_jet_lat_vals:
            lats_within_jet.append(current_lat)
        else:
            current_lat_is_within_range = False
        current_lat += amount_to_add_to_lat_val
    return lats_within_jet


def scale_lat_vals_with_quadratic_func(lats, speeds, scaled_lats):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2017) https://doi.org/10.1175/JCLI-D-16-0849.1

    Will downscale or upscale the resolution of latitude using a quadratic func

    Parameters
    ----------
    lats : xr.DataArray or array-like
        Array of latitude values
    speeds :  xr.DataArray or array-like
        Array of wind-speeds
    scaled_lats : array-like
        Array of scaled latitude values of a given resolution (see rescale_lat_resolution function)

    Returns
    ----------
    scaled_lat_vals : xr.DataArray or array-like
        Array of rescaled latitude values based scaled_lats
    """
    scaled_lat_vals = apply_quadratic_func(lats, speeds, scaled_lats)
    return scaled_lat_vals


def rescale_lat_resolution(lats, lat_resolution):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani 2017 https://doi.org/10.1175/JCLI-D-16-0849.1
    & Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1

    TODO: what if larger resolution

    Parameters
    ----------
    lats : xr.DataArray or array-like
        Array of latitude values
    lat_resolution : int or float
        Latitude resolution in degrees

    Returns
    ----------
    output : numpy.array
        Rescaled array of latitude values
    """
    return np.arange(min(lats), max(lats) + lat_resolution, lat_resolution)


def get_latitude_and_speed_where_max_ws_at_reduced_resolution(
    lats_and_ws, lat_resolution
):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2017) https://doi.org/10.1175/JCLI-D-16-0849.1

    Makes use of the quadratic func to scale latitude values

    Parameters
    ----------
    lats_and_ws : xr.DataArray or array-like
        Array of latitudes and windspeeds
    lat_resolution : int or float
        Latitude resolution in degrees

    Returns
    ----------
    output : numpy.array
        latitude and wind-speed value scaled by quadratic func
    """
    lats, ws = lats_and_ws
    #  Remove numpy.nan from list
    lats = [lat for lat in lats if not np.isnan(lat)]
    ws = [s for s in ws if not np.isnan(s)]
    #  Scale lats
    scaled_lats = rescale_lat_resolution(lats, lat_resolution)
    scaled_lat_vals = scale_lat_vals_with_quadratic_func(lats, ws, scaled_lats)
    decimal_places = general_utils.get_num_of_decimal_places(lat_resolution)
    max_speed_at_scaled_lat = np.max(scaled_lat_vals)
    return (
        round(scaled_lats[np.argmax(scaled_lat_vals)], decimal_places),
        max_speed_at_scaled_lat,
    )


def get_3_latitudes_and_speed_around_max_ws(row):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2017) https://doi.org/10.1175/JCLI-D-16-0849.1

    Will get the latitudes neighbouring to east, at and to west where the max windspeed is found

    Parameters
    --------------
    row : xr.DataArray or array-like
        Array containing latitude and wind-speed values

    Returns
    ----------
    neighbouring_lats : array-like
        array of 3 neighbouring latitude coordinates around where maximum windspeed is found in input (row)
    neighbouring_speeds : array-like
        array of 3 neighbouring winspeeds values around where maximum windspeed is found in input (row)
    """
    assert "lat" in row.coords, "'lat' needs to be in data.coords"

    lat_resolution = float(row["lat"][1] - row["lat"][0])
    lat_min, lat_max = float(row["lat"].min()), float(row["lat"].max())
    max_lat, _ = get_latitude_and_speed_where_max_ws(row)
    if np.isnan(max_lat):
        # occurs when no data in slice
        return (
            np.array([np.nan, np.nan, np.nan], dtype="float64"),
            np.array([np.nan, np.nan, np.nan], dtype="float64"),
        )
    neighbouring_lats = get_3_neighbouring_coord_values(
        max_lat, lat_resolution
    )
    neighbouring_lats = neighbouring_lats[
        (neighbouring_lats >= lat_min) & (neighbouring_lats <= lat_max)
    ]
    neighbouring_speeds = row.sel(lat=neighbouring_lats).data

    #  TODO: move work around to func
    if len(neighbouring_lats) < 3:
        if neighbouring_lats[0] == lat_min:
            neighbouring_lats = np.insert(neighbouring_lats, 0, np.nan)
            neighbouring_speeds = np.insert(neighbouring_speeds, 0, np.nan)
        if neighbouring_lats[-1] == lat_max:
            neighbouring_lats = np.append(neighbouring_lats, np.nan)
            neighbouring_speeds = np.append(neighbouring_speeds, np.nan)
    return (neighbouring_lats, neighbouring_speeds)


def get_3_neighbouring_coord_values(coord_val, coord_resolution):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2017) https://doi.org/10.1175/JCLI-D-16-0849.1

    TODO: add to JetStreamOccurenceAndCentreAlgorithm and ...

    Parameters
    ----------
    coord_val : int or float
        Central coord value to get neighbours from
    coord_resolution : int or float
        in degrees

    Returns
    ----------
    output : array-like
        array of 3 neighbouring latitude coordinates

    Usage
    ----------
    get_3_neighbouring_coord_values(45.0, 1.25)
    >>> np.array([43.75, 45.0, 46.25])
    """
    if not isinstance(coord_val, float) or not isinstance(
        coord_resolution, float
    ):
        coord_val = float(coord_val)
        coord_resolution = float(coord_resolution)

    return np.array(
        [coord_val - coord_resolution, coord_val, coord_val + coord_resolution]
    )


def quadratic_func(x, y):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2017) https://doi.org/10.1175/JCLI-D-16-0849.1


    Parameters
    ----------
    x : xr.DataArray or array-like
        Array 1
    y : xr.DataArray or array-like
        Array 2
    """
    p = np.polyfit(x, y, deg=2)
    return p


def apply_quadratic_func(x, y, vals):
    """
    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2017) https://doi.org/10.1175/JCLI-D-16-0849.1


    Parameters
    ----------
    x : xr.DataArray or array-like
        Array 1
    y : xr.DataArray or array-like
        Array 2
    vals : array-like
        Values to apply function to

    Returns
    ----------
    output : array-like
        quadratic function output
    """
    a, b, c = quadratic_func(x, y)
    return (a * vals**2) + (b * vals) + c


def make_empty_local_wind_maxima_data_var(data):
    """
    Will add a new data var of zeros for local wind maxima

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    TODO: add asserts
    """
    data["local_wind_maxima"] = (
        ("time", "plev", "lat", "lon"),
        np.zeros(
            (
                len(data["time"]),
                len(data["plev"]),
                len(data["lat"]),
                len(data["lon"]),
            )
        ),
    )
    return data


def get_empty_local_wind_maxima_data(data):
    """
    Will add a new data var of zeros for local wind maxima

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    data : xarray.Dataset
        Input data with (time, plev, lat, lon) dimensions

    Returns
    ----------
    data : xarray.Dataset
        Data containing zeros array of (time, plev, lat, lon) dimensions
    TODO: add asserts
    """
    data["local_wind_maxima"] = (
        ("time", "plev", "lat", "lon"),
        np.zeros(
            (
                len(data["time"]),
                len(data["plev"]),
                len(data["lat"]),
                len(data["lon"]),
            )
        ),
    )
    return data


def get_potential_local_wind_maximas_by_ws_threshold(
    ws_slice, ws_threshold=30
):
    """
    Will return a 2-d array of potential local windspeed maximas

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    ws_slice : xarray.Dataset
        Data slice of windspeed that has only lat and lon dims

    ws_threshold : int or float
        windspeed threshold to apply (default=30 ms-1)
    Returns
    ----------
    ws_slice : xarray.Dataset
        Data slice of windspeed (lat, lon only) with ws_threshold applied
    TODO: add checks
    """
    return ws_slice.where(lambda x: x > ws_threshold).fillna(0.0)


def get_local_wind_maxima_by_timeunit(row):
    """
    Get local wind maxima by timeunit (i.e. day)

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    row : xarray.Dataset
        Data of a single time unit containing windspeed (ws)

    Returns
    ----------
    row : xarray.Dataset
        Data of a single time unit containing 0 or 1 value (local wind maxima) for that time unit
    """
    if "local_wind_maxima" not in row.data_vars:
        raise ValueError("local_wind_maxima needs to be defined.")

    row = row.transpose("plev", "lat", ...)
    for lon in row["lon"]:
        current = row.sel(lon=lon)
        pot_local_maximas = get_potential_local_wind_maximas_by_ws_threshold(
            current["ws"], 30
        ).data
        ind_local_wind_maximas = general_utils.get_local_maxima(
            pot_local_maximas, axis=1
        )
        # Turn into 2-d numpy array
        ind_local_wind_maximas = np.array(
            [
                [arr1, arr2]
                for arr1, arr2 in zip(
                    ind_local_wind_maximas[0], ind_local_wind_maximas[1]
                )
            ]
        )
        for plev_ind, lat_ind in ind_local_wind_maximas:
            row["local_wind_maxima"].loc[
                dict(
                    plev=current["plev"].data[plev_ind],
                    lon=lon,
                    lat=current["lat"].data[lat_ind],
                )
            ] = 1.0
    return row


def get_number_of_timeunits_per_monthyear_with_local_wind_maxima(data):
    """
    Will resample by each month and return number of timeunits (i.e. day) with local wind maxima

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    data : xarray.Dataset
        Input data with local wind maxima values by time unit (i.e. day)

    Returns
    ----------
    data : xarray.Dataset
        Data containing zeros array of (time, plev, lat, lon) dimensions

    """
    data = (
        data["local_wind_maxima"]
        .resample(time="MS")
        .sum()
        .rename({"time": "monthyear"})
    )
    return data


def subdivide_local_wind_maxima_into_stj_pfj(data):
    """
    Subdivide the local_wind_maxima values into the Subtropical Jet (STJ) and Polar Front Jet (PFJ) based on pg. 2709

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    data : xarray.Dataset
        Input data with local wind maxima values

    Returns
    ----------
    data : xarray.Dataset
        data with polar_front_jet and subtropical_jet subdivisions
    """
    DJF_STJ = data.sel(
        monthyear=data.monthyear.dt.month.isin([1, 2, 12]), lat=slice(15, 40)
    )["local_wind_maxima"]
    MAM_SON_PFJ = data.sel(
        monthyear=data.monthyear.dt.month.isin([3, 4, 5, 9, 10, 11]),
        lat=slice(10, 70),
    )["local_wind_maxima"]
    JJA_PFJ = data.sel(
        monthyear=data.monthyear.dt.month.isin([6, 7, 8]), lat=slice(30, 60)
    )["local_wind_maxima"]
    SH_STJ = data.sel(lat=slice(-40, -15))["local_wind_maxima"]
    SH_PFJ = data.sel(lat=slice(-70, -41))["local_wind_maxima"]
    data["polar_front_jet"] = xr.merge([MAM_SON_PFJ, JJA_PFJ, SH_PFJ])[
        "local_wind_maxima"
    ]
    data["subtropical_jet"] = xr.merge([DJF_STJ, SH_STJ])["local_wind_maxima"]
    return data


def run_jet_occurence_and_centre_alg_on_one_day(row, occurence_ws_threshold):
    """
    Component of method from Kuang et al (2014) https://doi.org/10.1007/s00704-013-0994-x

    Runs JetStreamCoreIdentificationAlgorithm method on a single day

    Parameters
    ----------
    row : xarray.Dataset
        Single time-unit of dataset containing u- and v-component wind
    occurence_ws_threshold : int or float
        Threshold used to identify a jet-stream occurence point

    Returns
    ----------
    occ_alg.output_data : xarray.Dataset
        Data with jet occurence and centre points (1 for occurence, 2 for centre)
    """
    occ_alg = JetStreamOccurenceAndCentreAlgorithm(row, occurence_ws_threshold)
    occ_alg.run()
    return occ_alg.output_data


class JetStreamOccurenceAndCentreAlgorithm:
    """
    Component of method from Kuang et al (2014) https://doi.org/10.1007/s00704-013-0994-x
    """

    def __init__(self, data, occurence_ws_threshold=30):
        """
        Parameters
        ----------
        data : xarray.Dataset
            Single time-unit of dataset containing u- and v-component wind. Needs to be 2-dimensions (lat, lon)
        occurence_ws_threshold : int or float
            Threshold used to identify a jet-stream occurence point
        """
        try:
            assert occurence_ws_threshold > 0
        except Exception as e:
            raise ValueError(
                "Occurence wind-speed threshold needs to be more than 0"
            ) from e

        # Load in data as a pressure level 2d wind-speed slice
        self.plev_ws_slice = windspeed_utils.PressureLevelWindSpeedSlice(
            data
        ).values
        self.plev_ws_slice["jet_ocurrence1_jet_centre2"] = self.plev_ws_slice[
            "ws"
        ].copy()
        self.plev_ws_slice["jet_ocurrence1_jet_centre2"] = self.plev_ws_slice[
            "jet_ocurrence1_jet_centre2"
        ].where(lambda x: x >= occurence_ws_threshold)
        self._jet_occurence = self.plev_ws_slice
        self._lat_resolution = float(
            self.plev_ws_slice["lat"][1] - self.plev_ws_slice["lat"][0]
        )
        self._lon_resolution = float(
            self.plev_ws_slice["lon"][1] - self.plev_ws_slice["lon"][0]
        )

        # make_duplicate data for output
        self.output_data = self._jet_occurence.copy(deep=True)

        # Initialise lists needed for search algorithm
        self._all_coords = []
        self._lats_with_3 = []
        self._lats_for_search = []
        self._jet_centres = []

        # Needed to keep track of algorithm
        self.algorithm_has_run = False

    @classmethod
    def run_algorithm(cls, data):
        return cls(data).run()

    def run(self):
        self._get_all_coords_of_jet_occurence()
        self._all_coords_arr = np.array(self._all_coords)
        # Get a counter of all the latitude coordinates
        # TODO: add tests to see if this works
        try:
            self._count_lats = collections.Counter(self._all_coords_arr[:, 0])
        except Exception as e:
            print(e)
            self._count_lats = {}
        self._get_all_lats_of_jet_centre_for_search()
        self._calc_jet_centre_points()
        self._get_jet_centre_data()
        self._label_jet_occurence()

        self.algorithm_has_run = True

    def _get_jet_centre_data(self):
        """
        Calculates jet-stream centres based on if one jet-stream occurence grid
        is surrounded by 8 cells of jet-stream occurence (default is 30 m/s)
        """
        # TODO: there's got to be a quicker way
        for centre in self._jet_centres:
            self.output_data["jet_ocurrence1_jet_centre2"].loc[
                dict(lat=centre[0], lon=centre[1])
            ] = 2

    def _get_all_coords_of_jet_occurence(self):
        for val in self._jet_occurence["jet_ocurrence1_jet_centre2"].notnull():
            if val.any():
                for sub_val in val:
                    if sub_val:
                        self._all_coords.append(
                            [float(sub_val["lat"]), float(sub_val["lon"])]
                        )

    def _get_all_lats_of_jet_centre_for_search(self):
        """
        Will get all latitudes that 'may' be a jet-stream centre-point
        i.e. where latitude appears at least three times for 3*3 lat/lon grid.
        NOTE: speeds up calculation as less values are searched through
        """
        # Step 1. look for all latitudes with at least 3 occurences - 3*3 grid
        self._get_all_latitudes_that_occur_at_least_three_times()
        # Step 2. Check if the latitudes above and below the lats with 3 values
        # are present i.e. Y component of for 3*3
        self._get_all_latitudes_available_in_3by3_grid()

    def _get_all_latitudes_that_occur_at_least_three_times(self):
        for lat in self._count_lats.items():
            if lat[1] >= 3:
                self._lats_with_3.append(lat[0])

    def _get_all_latitudes_available_in_3by3_grid(self):
        for lat in self._lats_with_3:
            if (
                lat - self._lat_resolution in self._lats_with_3
                and lat + self._lat_resolution in self._lats_with_3
            ):
                self._lats_for_search.append(lat)

    def _calc_jet_centre_points(self):
        """
        Will return a list of the coordinates for all jet-centre points
        """
        for lat in self._lats_for_search:
            coords_to_search = self._all_coords_arr[
                np.where(self._all_coords_arr[::, 0] == lat)
            ]
            for coord in coords_to_search:
                if (
                    coord[0] == 0
                    or coord[1] == 0
                    and 360 - self._lon_resolution
                    not in coords_to_search[::, 1]
                ):
                    continue
                # check if coord is jet centre point i.e. 9*9 all above 30
                lat_grid_vals = np.arange(
                    coord[0] - self._lat_resolution,
                    coord[0] + self._lat_resolution + 0.1,
                    self._lat_resolution,
                )
                lon_grid_vals = np.arange(
                    coord[1] - self._lon_resolution,
                    coord[1] + self._lon_resolution + 0.1,
                    self._lon_resolution,
                )
                matrix_vals_to_check = np.array(
                    np.meshgrid(lat_grid_vals, lon_grid_vals)
                ).T.reshape(-1, 2)
                matrix_vals_to_check = (
                    matrix_vals_to_check % 360
                )  # loop around
                add_coord = True
                for val in matrix_vals_to_check:
                    if not val.tolist() in self._all_coords:
                        add_coord = False
                        break
                if add_coord:
                    self._jet_centres.append(coord)

    def _label_jet_occurence(self):
        """
        Will label all non 2 values of windspeed for the occurence
        Used in Kuang et al. 2014
        """
        self.output_data["jet_ocurrence1_jet_centre2"] = self.output_data[
            "jet_ocurrence1_jet_centre2"
        ].where(lambda x: np.isfinite(x), 0)
        self.output_data["jet_ocurrence1_jet_centre2"] = self.output_data[
            "jet_ocurrence1_jet_centre2"
        ].where(lambda x: ((x == 0) | (x == 2)), 1)


def get_jet_lat_and_speed_using_parabola_by_day(data_row):
    """
    Will get jet latitude and speed by fitting a parabola and taking maximum

    Component of method from Barnes & Polvani (2015) http://journals.ametsoc.org/doi/10.1175/JCLI-D-14-00589.1

    Parameters
    ----------
    data_row : xarray.DataArray
        Data of single time unit containing u-component wind

    Returns
    ----------
    data_row : xarray.DataArray
        Data of single time unit containing jet_lat and jet_speed variables
    """
    fitted_parabola = fit_parabola(data_row["lat"].data, data_row["ua"].data)
    ind_of_max = fitted_parabola.argmax()
    data_row["jet_lat"] = float(data_row.lat[ind_of_max])
    data_row["jet_speed"] = fitted_parabola.max()
    return data_row


def fit_parabola(x, y):
    """
    Fits a parabola
    TODO: check if correct

    Component of method from Barnes & Polvani (2015) http://journals.ametsoc.org/doi/10.1175/JCLI-D-14-00589.1
    """
    coeff, cov = scipy.optimize.curve_fit(parabola, x, y)
    fitted_parabola = parabola(x, coeff[0], coeff[1], coeff[2])
    return fitted_parabola


def parabola(x, a, b, c):
    """
    Parabola
    Component of method from Barnes & Polvani (2015) http://journals.ametsoc.org/doi/10.1175/JCLI-D-14-00589.1
    """
    return a * x**2 + b * x + c


def calc_meridional_circulation_index(data):
    """
    Component of method from Francis and Vavrus (2015) https://doi.org/10.1088/1748-9326/10/1/014005
    Calculates the Meridional Circulation Index (MCI)
    When MCI = 0, the wind is purely zonal, and when MCI= 1 (-1), the flow is
    from the South (North).

           v * abs(v)
    MCI =  ――――――――――
           u**2 * v**2

    NOTE: The paper is not clear about whether the absolute value for MCI
    is taken instead thus 0-1

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component and v-component wind-speed values

    Returns
    ----------
    output : xarray.DataArray
        Array of Meridional Circulation Index (MCI) values
    """
    assert (
        "ua" in data.variables and "va" in data.variables
    ), "Cannot compute metric. 'ua' and/or 'va' not found in data"
    return data["va"] * abs(data["va"]) / (data["ua"] ** 2 + data["va"] ** 2)


def get_latitude_circle_linestring(latitude, lon_min, lon_max):
    """
    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309
    Will return a linestring of a latitude circle

    Parameters
    ----------
    latitude : int or float
        given latitude to calculate circle from
    lon_min : int or float
        Minimum longitude for circle to extend to
    lon_max : int or float
        Maximum longitude for circle to extend to

    Returns
    ----------
    circle : shapely.geometry.LineString
        Linestring of latitude circle around a hemisphere
    """
    vals = np.column_stack(
        (
            np.arange(lon_min, lon_max + 0.1, 0.5),
            np.array([latitude] * len(np.arange(lon_min, lon_max + 0.1, 0.5))),
        )
    )
    circle = shapely.geometry.LineString(vals)
    return circle


def get_sinuosity_of_zonal_mean_zg(row, latitude_circle):
    """
    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309
    Works on a grouped data set and will calculate sinuosity of zonal mean
    geopotential (ZG) contour compared to a latitude circle

    Parameters
    ----------
    row : xarray.Dataset
        Data containing geopotential height (zg) and zonal_mean_zg_30Nto70N values
    latitude_circle : shapely.geometry.LineString
        Linestring of latitude circle around a hemisphere

    Returns
    ----------
    row : xarray.Dataset
        Data containing sinuousity value (determined by calc_great_circle_sinuosity function)
    """
    row["sinuosity"] = calc_great_circle_sinuosity(
        get_one_contour_linestring(
            row["zg"], row["zonal_mean_zg_30Nto70N"].data
        ),
        latitude_circle,
    )
    return row


def get_one_contour_linestring(dataarray, contour_level):
    """
    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Returns a linestring or multi-linestring of a given contour

    Parameters
    ----------
    dataarray : xarray.DataArray
        Array of Geopotential height (zg) values to calculate contour from
    contour_level :
        Value with which to calculate a contour from geopotential height

    Returns
    ----------
    contour_line : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Contour line of geopotential height (zg) a given contour
    """
    assert isinstance(
        dataarray, xr.DataArray
    ), "Data needs to be type xr.DataArray"
    assert (
        "lat" in dataarray.coords and "lon" in dataarray.coords
    ), "Data array needs to have latitude and longitude coords"
    one_contour = dataarray.plot.contour(levels=[contour_level])
    matplotlib.pyplot.close()
    if len(one_contour.allsegs[0]) > 1:
        try:
            contour_line = shapely.geometry.MultiLineString(
                (one_contour.allsegs[0])
            )
        except ValueError as ve:
            print(ve)
            return np.nan
    else:
        contour_line = shapely.geometry.LineString((one_contour.allsegs[0][0]))
    return contour_line


def calc_total_great_circle_distance_along_line(line):
    """
    Returns the total great circle (haversine) distance along a linestring
    or multilinestring

    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Parameters
    ----------
    line : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Line to calculate great circle distance from

    Returns
    ----------
    total_distance : int or float
        Total distance in degrees or m? TODO

    """
    total_distance = 0
    if isinstance(line, shapely.geometry.multilinestring.MultiLineString):
        for i, _ in enumerate(line):
            total_distance += (
                general_utils.get_great_circle_distance_along_linestring(
                    shapely.geometry.LineString((line[i]))
                )
            )
    elif isinstance(line, shapely.geometry.LineString):
        total_distance += (
            general_utils.get_great_circle_distance_along_linestring(line)
        )
    else:
        return np.nan
    return total_distance


def calc_great_circle_sinuosity(line1, line2):
    """
    Calculates sinuosity by comparing the great circle distance between
    two (multi-)linestrings

    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Parameters
    ----------
    line1 : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Line to calculate great circle sinuosity from

    line2 : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Line to calculate great circle distance from

    Returns
    ----------
    sinuosity : int or float
        Sinuosity value between length of two lines
    """
    return calc_total_great_circle_distance_along_line(
        line1
    ) / calc_total_great_circle_distance_along_line(line2)


def assign_jet_lat_speed_to_ten_day_mean_data(ten_day_mean):
    """
    Joins the values for latitude and windspeed of the point of maximum windspeed for a given sector to the ten_day_mean data

    Component of method from Barnes & Simpson 2017 https://doi.org/10.1175/JCLI-D-17-0299.1

    Parameters
    ----------
    ten_day_mean : xarray.Dataset
        Data containing u-component wind and resampled by 10 days

    Returns
    ----------
    ten_day_mean : xarray.Dataset
        Data with jet latitude and windspeed
    """
    max_lats_ws = np.array(
        list(map(get_latitude_and_speed_where_max_ws, ten_day_mean["ua"]))
    )
    ten_day_mean = assign_jet_lat_and_speed_to_data(ten_day_mean, max_lats_ws)
    return ten_day_mean


def assign_ten_day_average_jet_lat_speed_to_data(data, ten_day_mean):
    """
    Joins the values for latitude and windspeed of the point of maximum windspeed to the main data

    Component of method from Barnes & Simpson 2017 https://doi.org/10.1175/JCLI-D-17-0299.1

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component wind

    Returns
    ----------
    output : xarray.Dataset
        Data with jet latitude and windspeed
    """
    return data.assign(
        {
            "jet_lat": (("10_day_average"), ten_day_mean["jet_lat"].data),
            "jet_speed": (("10_day_average"), ten_day_mean["jet_speed"].data),
        }
    )


def cubic_spline_interpolation(x, y):
    """
    Component of method from Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1
    Runs a cubic spline interpolation using 2 equal-length arrays

    Parameters
    ----------
    x : xarray.DataArray or array-like
        array 1 to run cubic spline interpolation
    y : xarray.DataArray or array-like
        array 2 to run cubic spline interpolation

    Returns
    ----------
    cubic_interpolation : np.array
        Cubic spline interpolation
    """
    return scipy.interpolate.interp1d(
        x, y, kind="cubic", fill_value="extrapolate"
    )


def run_cubic_spline_interpolation_to_get_max_lat_and_ws(
    data, lat_resolution, ua_col="ua"
):
    """
    Component of method from Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1
    Runs a cubic spline interpolation to find maximum latitude and maximum windspeed at a given resolution of latitude

    Parameters
    ----------
    data : xr.Dataset
        must contain coords lat
    lat_resolution : int or float
        Latitude resolution in degrees for cubic spline interpolation
    ua_col : str
        u-component windspeed column (default='ua')

    Returns
    ----------
    max_lat : float
        latitude with the maximum wind-speed for time period as deterimined by cubic spline interpolation to a given lat resolution
    max_ws : float
        maximum wind-speed at the given latitude for time period
    """
    scaled_lats = rescale_lat_resolution(data["lat"], lat_resolution)
    csi = cubic_spline_interpolation(data["lat"], data[ua_col])
    interpolated_ws = csi(scaled_lats)
    max_lat = scaled_lats[np.argmax(interpolated_ws)]
    max_ws = max(interpolated_ws)
    return max_lat, max_ws


def run_cubic_spline_interpolation_for_each_unit_of_climatology_to_get_max_lat_and_ws(
    data, lat_resolution, time_col
):
    """
    Component of method from Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1
    Runs a cubic spline interpolation to find maximum latitude and maximum windspeed at a given resolution of latitude for each

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component wind-speed. In the case of Bracegirdle et al. 2018, uses a seasonal and annual mean/climatology
    resolution : int or float
        Latitude resolution in degrees for cubic spline interpolation
    time_col : str
        Column name containing period (i.e. season, year, month, etc.)

    Returns
    ----------
    max_lats : list
        list of latitudes with the maximum wind-speed
    max_ws : list
        list of maximum wind-speed at the given latitudes

    Usage
    ----------
    seasonal_climatology = general_utils.get_climatology(data, "season")
    seasonal_zonal_mean = seasonal_climatology.mean("lon")
    (
        seasonal_max_lats,
        seasonal_max_ws,
    ) = jetstream_metrics_utils.run_cubic_spline_interpolation_for_each_unit_of_climatology_to_get_max_lat_and_ws(
        seasonal_zonal_mean, resolution=0.075, time_col="season"
    )
    """
    max_lats = []
    max_ws = []
    for period in data[time_col].data:
        period_data = data.sel({time_col: period})
        lat, ws = run_cubic_spline_interpolation_to_get_max_lat_and_ws(
            period_data, lat_resolution=lat_resolution
        )
        max_lats.append(lat)
        max_ws.append(ws)
    return max_lats, max_ws


def calc_centroid_jet_lat_from_zonal_mean(zonal_mean, area_by_lat):
    """
    Component of method from Ceppi et al (2018) https://doi.org/10.1175/JCLI-D-17-0323.1

    Will get the centroid latitude of the u-component wind using:

                integral(60deg, 30deg)(zonal_u**2*lat) dlat
    jet_lat =   ------------------------------
                integral(60deg, 30deg)(zonal_u**2) dlat

    Parameters
    ----------
    zonal_mean : xarray.Dataset
        zonally-average data containing u-component wind data
    area_by_lat : xarray.DataArray
        Information on area in each latitude (i.e. m2 per latitude)
    """
    u_hat_by_lat = zonal_mean["ua"] ** 2 * zonal_mean["lat"] * area_by_lat
    u_hat = zonal_mean["ua"] ** 2 * area_by_lat
    return u_hat_by_lat.sum("lat") / u_hat.sum("lat")


def get_moving_averaged_smoothed_jet_lats_for_one_day(data_row):
    """
    Component of method from Kerr et al. (2020) https://onlinelibrary.wiley.com/doi/10.1029/2020JD032735
    Parameters
    ----------
    data_row : xarray.DataArray
        Data containing u-component windspeed of one time unit

    Returns
    ----------
    smoothed_jet_lat : xarray.Dataset
        Data containing jet-stream position
    """
    data_row["jet_lat_by_lon"] = get_jet_lat_by_lon(data_row["ua"])
    data_row[
        "smoothed_jet_lats"
    ] = smooth_jet_lat_across_lon_with_rectangular_pulse(
        data_row["jet_lat_by_lon"], width_of_pulse=10
    )
    return data_row


def smooth_jet_lat_across_lon_with_rectangular_pulse(
    jet_lat_data, width_of_pulse
):
    """
    Smooth jet position (jet latitude) by carrying out a convolution with a rectangular pulse

    Component of method from Kerr et al. (2020) https://onlinelibrary.wiley.com/doi/10.1029/2020JD032735

    Parameters
    ----------
    jet_lat_data : xarray.DataArray
        Data detailing the latitude of jet-stream at each longitude of one time unit
    width_of_pulse : float or int
        Width to make rectangular pulse (likely to be in  'degrees_east')

    Returns
    ----------
    filtered_data : xarray.DataArray
        Data detailing the latitude of jet-stream that has been smoothed using rectangular pulse

    """
    sig = jet_lat_data["lon"] % jet_lat_data < width_of_pulse
    filtered_data = jet_lat_data[sig]
    return filtered_data


def get_jet_lat_by_lon(data_row):
    """
    Gets all the latitudes of maximum wind-speeds by longitude

    Component of method from Kerr et al. (2020) https://onlinelibrary.wiley.com/doi/10.1029/2020JD032735

    Parameters
    ----------
    data_row : xarray.DataArray
        Data containing u-component windspeed of one time unit

    Returns
    ----------
    output : xarray.DataArray
        Data detailing the latitude of jet-stream at each longitude of one time unit
    """
    return data_row["lat"][data_row.argmax("lat")]
