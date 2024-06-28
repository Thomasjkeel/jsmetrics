# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to
    identify or classify jet-stream in the literature.

    This file is in order of publish year of the metrics
    (see details_for_all_metrics.py)

    Functions and Classes ordered by the year of the paper that uses each metric component first
"""

# imports
import warnings
import numpy as np
import xarray as xr
import scipy.fftpack
import scipy.interpolate
import scipy.optimize
from jsmetrics.utils import data_utils

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def calc_atmospheric_mass_at_kPa(pressure, gravity=9.81, atmospheric_area=5.1e8):
    """
    Will calculate the atmospheric mass at a given pressure level.

    Radius of earth (R) = 6.372E3 km;
    Surface area of earth = 4 Pi R^2 = 5.1E8 km^2

    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    Get atmospheric mass at one pressure level.

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


def get_weighted_average_at_one_pressure_level(data, pressure_level, atm_mass, ws_col):
    """
    Get weighted average wind speed at one pressure level.

    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Parameters
    ----------
    data : xarray.Dataset
        Data containing windspeed ('ws') or u-, v-component wind
    pressure_level : int or float
        Pressure level to subset
    atm_mass : int or float
        Atmospheric mass at given hPa pressure level
    ws_col : string
        Name of column to calculate weighted average from (e.g. 'ws', 'ua', 'va')

    Returns
    ----------
    output : xarray.Dataset
        Data with weighted average at a single pressure level
    """
    return atm_mass * (data[ws_col].sel(plev=pressure_level))


def get_mass_weighted_average_wind(data, ws_col, plev_flux=False):
    """
    Get mass-weighted average wind-speed from 'u', 'v' component wind or 'wind vector'.

    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    data_utils.check_at_least_n_plevs_in_data(data, n_plevs=2)
    sum_weighted_ws = None
    (
        hPa_multiplying_factor,
        _,
    ) = get_hPa_or_mbar_and_Pa_or_bar_multiplying_factors_for_plev(data)
    for plev_step in data["plev"].data:
        plev_hPa = plev_step * hPa_multiplying_factor
        atm_mass = get_atm_mass_at_one_hPa(plev_hPa)
        weighted_average = get_weighted_average_at_one_pressure_level(
            data, plev_step, atm_mass, ws_col
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
    Sum atmospheric mass for given pressure levels.

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
    data_utils.check_at_least_n_plevs_in_data(data, n_plevs=2)
    (
        hPa_multiplying_factor,
        _,
    ) = get_hPa_or_mbar_and_Pa_or_bar_multiplying_factors_for_plev(data)
    for plev_level in data["plev"].data:
        plev_hPa = plev_level * hPa_multiplying_factor
        atm_mass = get_atm_mass_at_one_hPa(plev_hPa)
        sum_atm_mass += atm_mass
    return sum_atm_mass


def get_hPa_or_mbar_and_Pa_or_bar_multiplying_factors_for_plev(data):
    """
    Will get the multiplying factor for plev units.
    E.g. if data contains plev with the unit: Pa the two multiplying factors will be 0.01 and 1
    E.g. if data contains plev with the unit: hPa the two multiplying factors will be 1 and 100

    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Parameters
    ----------
    data : xr.DataArray or xr.Dataset
        Data containing plev coord and units for plev

    Returns
    ----------
    hPa_multiplying_factor : int
        Factor to multiply plev units to convert them to hPa
    Pa_multiplying_factor : int
        Factor to multiply plev units to convert them to Pa

    """
    hPa_multiplying_factor = 1
    Pa_multiplying_factor = 1
    plev_units = data_utils.check_plev_units(
        data, ["Pa", "hPa", "mbar", "millibars", "bars"]
    )

    if plev_units in ["hPa", "mbar", "millibars"]:
        Pa_multiplying_factor *= 100
    else:
        hPa_multiplying_factor /= 100

    return hPa_multiplying_factor, Pa_multiplying_factor


def calc_mass_weighted_average(data, ws_col):
    """
    Calculate mass weighted average wind-speed.
    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    Calculate mass-flux weighted pressure.

    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614

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
    Calculate mass-flux weighted latitude.

    Component of method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614
    NOTE: Problem with including 1000 hPa

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
        lat_data = sub_data.sel(lat=lat, method="nearest")
        lat_sum_weighted_ws = get_mass_weighted_average_wind(lat_data, ws_col)
        if sum_weighted_lat_flux is None:
            sum_weighted_ws_by_lat = lat_sum_weighted_ws
            sum_weighted_lat_flux = lat_sum_weighted_ws * lat
        else:
            sum_weighted_ws_by_lat += lat_sum_weighted_ws
            sum_weighted_lat_flux += lat_sum_weighted_ws * lat
    mass_flux_weighted_latitude = sum_weighted_lat_flux / sum_weighted_ws_by_lat
    return mass_flux_weighted_latitude


def get_climatology(data, freq):
    """
    Makes a climatology at given interval (i.e. days, months, season)

    Parameters
    ----------
    data : xarray.Dataset
        data with regular time stamp
    freq : str
        'day', 'month' or 'season'

    Returns
    ----------
    climatology : xarray.Dataset
        Climatology of a given frequency

    Usage
    ----------
    climatology = get_climatology(data, 'month')

    """
    climatology = data.groupby("time.%s" % (freq)).mean("time")
    return climatology


def calc_low_pass_weights(window, cutoff):
    """
    Calculate weights for a low pass Lanczos filter.
    A low-pass filter removes short-term random fluctations in a time series

    Component of method from Woollings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    Will carry out Lanczos low-pass filter

    Component of method from Woollings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    if window_size <= filter_freq or filter_freq <= 0 or window_size <= 0:
        raise ValueError(
            "Lanczos filter frequency cannot be larger than window size and both need to be more than 0"
        )

    if not xr.infer_freq(dataarray["time"]) == "D":
        print(
            f"Warning: The 'apply_lanczos_filter' function was built to work on a datetime index of freq 'D'. Frequency in data is '{xr.infer_freq(dataarray['time'])}'"
        )

    assert isinstance(dataarray, xr.DataArray), "Input data needs to be a data array"

    try:
        start_date = dataarray["time"].astype(np.datetime64)[0]
        filter_end_date = start_date + np.timedelta64(
            filter_freq, "D"
        )  # add filter end date for check if data is outside the filter size
        window_end_date = start_date + np.timedelta64(
            window_size, "D"
        )  # add window end date for check if data is outside the filter size
        end_date = dataarray["time"].astype(np.datetime64)[-1]
    except (ValueError, TypeError):
        # makes assumption that if times not coerced into numpy datetime then is 360
        start_date = dataarray["time"].values[0]
        end_date = dataarray["time"].values[-1]
        if not hasattr(start_date, "calendar"):
            raise ValueError(
                f"Error: datetime type inputted ({start_date}) cannot be coerced into numpy.datetime64 and is not in cftime format."
            )
        if start_date.calendar == "360_day":
            datetime_added_func_to_use = data_utils.add_num_of_days_to_360Datetime
        elif start_date.calendar == "noleap":
            datetime_added_func_to_use = data_utils.add_num_of_days_to_NoLeapDatetime
        else:
            raise ValueError(
                f"Error: datetime type inputted ({start_date}) cannot be coerced into numpy.datetime64 and is not in cftime format."
            )
        filter_end_date = datetime_added_func_to_use(
            start_date, num_of_days_to_add=filter_freq
        )
        window_end_date = datetime_added_func_to_use(
            start_date, num_of_days_to_add=window_size
        )

    # add filter day for check if data is outside the filter size
    if window_end_date > end_date or filter_end_date > end_date:
        raise ValueError(
            f"Time series is too short to apply {window_size} day window for Lanczos filter frequency of {filter_freq} days"
        )

    #  As the data may not actually be in day format, it will return the number of values within window in data i.e. if in monthly format but filter window is 60, then actual filter size will be 2
    actual_filter_freq = dataarray.sel(
        time=slice(start_date, filter_end_date)
    ).time.size
    actual_window_size = dataarray.sel(
        time=slice(start_date, window_end_date)
    ).time.size
    lanczos_weights = calc_low_pass_weights(actual_window_size, 1 / actual_filter_freq)
    lanczos_weights_arr = xr.DataArray(lanczos_weights, dims=["window"])
    window_cons = (
        dataarray.rolling(time=len(lanczos_weights_arr), center=True)
        .construct("window")
        .dot(lanczos_weights_arr)
    )
    return window_cons


def get_latitude_and_speed_where_max_ws(data_row):
    """
    Will return the latitude and windspeed at the index of maximum wind speed
    from a row of data
    TODO: expects only one cell to be max

    Component of method from Woollings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    Barnes & Simpson (2017) https://doi.org/10.1175/JCLI-D-17-0299.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687

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
            print(
                "'get_latitude_and_speed_where_max_ws method': Warning: more than one max value found, picking the first occurence!"
            )
            max_ws = max_ws[0].max()
        lat_at_max = float(max_ws["lat"])
        speed_at_max = float(max_ws.data)
        return lat_at_max, speed_at_max
    else:
        return np.nan, np.nan


def calc_latitude_and_speed_where_max_ws(data, var_col="ua"):
    """
    Will return the latitude and windspeed at the index of maximum wind speed
    from a row of data
    TODO: expects only one cell to be max

    Component of method from Woollings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    Barnes & Simpson (2017) https://doi.org/10.1175/JCLI-D-17-0299.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687

    Parameters
    ----------
    dataarray : xarray.dataset
        Dataset containing u-component wind

    Returns
    ----------
    jet_lats :
        Latitude of maximum windspeed
    jet_speeds : int or float
        Speed at latitude of maximum windspeed
    """
    inds = np.nanargmax(data[var_col], axis=data[var_col].get_axis_num("lat"))
    jet_lats = data[var_col]["lat"][inds]
    data["jet_lat"] = (("time",), jet_lats.data)
    data["jet_speed"] = data[var_col].max(dim="lat")
    # cleans up missing values
    data = data.dropna("time", subset=["jet_lat", "jet_speed"])
    return data


def assign_jet_lat_and_speed_to_data(
    data, max_lat_ws, max_lats_col="jet_lat", max_ws_col="jet_speed"
):
    """
    Will return a data array with the maximum windspeed and latitude of that maximum wind speed

    Component of method from Woollings et al (2010) http://dx.doi.org/10.1002/qj.625
    & Barnes & Simpson 2017 https://doi.org/10.1175/JCLI-D-17-0299.1

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
    data_with_max_lats_ws[max_lats_col] = data_with_max_lats_ws[max_lats_col].astype(
        float
    )
    data_with_max_lats_ws[max_ws_col] = data_with_max_lats_ws[max_ws_col].astype(float)
    return data_with_max_lats_ws


def apply_low_freq_fourier_filter(data, highest_freq_to_keep):
    """
    Carries out a Fourier transform for filtering keeping only low frequencies.

    Component of method from Woollings et al (2010) http://dx.doi.org/10.1002/qj.625

    ADAPTED FROM: https://scipy-lectures.org/intro/scipy/auto_examples/plot_fftpack.html

    Parameters
    ----------
    data : np.array - 1-d
        time series data at regular intervals
    highest_freq_to_keep : int
        highest frequency to keep in the fourier transform expression
        NOTE: starts at 0, so highest_freq_to_keep=1 will only keep the constant and first expresion

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


def assign_filtered_lats_and_ws_to_data(data, filtered_max_lats, filtered_max_ws, dim):
    """
    Assigns the filtered data back to the returned dataset

    Component of method from Woollings et al (2010) http://dx.doi.org/10.1002/qj.625

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


def calc_jet_width_for_one_day(data_row, jet_lat, jet_speed):
    """
    Calculates jet width on a given data array using the calculated jet lat and speed at that lat

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    possible_surrounding_jet_vals = get_possible_surrounding_jet_lat_vals_for_one_day(
        data_row, jet_speed
    )
    if possible_surrounding_jet_vals.size == data_row["lat"].size:
        print("No jet-width determined")
        return np.nan
    possible_surrounding_jet_lat_vals = possible_surrounding_jet_vals["lat"]
    jet_width = get_width_of_jet_from_surrounding_lat_vals(
        jet_lat, possible_surrounding_jet_lat_vals, lat_resolution
    )
    return jet_width


def get_possible_surrounding_jet_lat_vals_for_one_day(one_day_wind_data, jet_speed):
    """
    Returns array of windspeed values that are within half speed of the jet_speed given.

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    Get latitude width of jet around the given jet latitude value. Width given in units of latitude

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    Returns all the latitudes in a continuous group (i.e. within one unit of lat resolution) around selected jet (jet_lat).

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    Creates a list of latitudes that are within a given distance from the jet_lat and only in one direction i.e. increasing or decreasing

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

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
    Will downscale or upscale the resolution of latitude using a quadratic func.

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687

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


def get_latitude_and_speed_where_max_ws_at_reduced_resolution(
    lats_and_ws, lat_resolution
):
    """
    Get latitude and speed where max windspeed at a reduced resolution.
    Makes use of the quadratic func to scale latitude values.

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687

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
    scaled_lats = data_utils.rescale_lat_resolution(lats, lat_resolution)
    scaled_lat_vals = scale_lat_vals_with_quadratic_func(lats, ws, scaled_lats)
    decimal_places = data_utils.get_num_of_decimal_places(lat_resolution)
    max_speed_at_scaled_lat = np.max(scaled_lat_vals)
    return (
        round(scaled_lats[np.argmax(scaled_lat_vals)], decimal_places),
        max_speed_at_scaled_lat,
    )


def get_3_latitudes_and_speed_around_max_ws(row):
    """
    Will get the latitudes neighbouring to east, at and to west where the max windspeed is found

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687

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

    max_lat, _ = get_latitude_and_speed_where_max_ws(row)
    if np.isnan(max_lat):
        # occurs when no data in slice
        return (
            np.array([np.nan, np.nan, np.nan], dtype="float64"),
            np.array([np.nan, np.nan, np.nan], dtype="float64"),
        )
    max_lat_argwhere = int(np.argwhere(row["lat"].data == max_lat))

    if max_lat_argwhere == 0:
        neighbouring_lats_ind = [max_lat_argwhere, max_lat_argwhere + 1]
    elif max_lat_argwhere == row["lat"].data.size - 1:
        neighbouring_lats_ind = [max_lat_argwhere - 1, max_lat_argwhere]
    else:
        neighbouring_lats_ind = [
            max_lat_argwhere - 1,
            max_lat_argwhere,
            max_lat_argwhere + 1,
        ]
    neighbouring_lats = row["lat"].isel(lat=neighbouring_lats_ind).data
    neighbouring_speeds = row.isel(lat=neighbouring_lats_ind).data

    #  add nan value to edges
    if len(neighbouring_lats) < 3:
        neighbouring_lats = add_nan_value_to_arr_if_not_len_enough(
            neighbouring_lats, len_required=3
        )
        neighbouring_speeds = add_nan_value_to_arr_if_not_len_enough(
            neighbouring_speeds, len_required=3
        )
    return (neighbouring_lats, neighbouring_speeds)


def add_nan_value_to_arr_if_not_len_enough(arr, len_required):
    """
    Add NaN value to array if any values in the array so it is a certain length

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687

    Parameters
    ----------
    arr : array-like
        Array to check.
    len_required : int
        length of array required

    Returns
    ----------
    arr : array-like
        new array

    """
    assert (
        len(arr) <= len_required
    ), "length of input array needs to be less or same as len_required"
    if len(arr) != len_required:
        arr = np.append(arr, np.nan)
        return add_nan_value_to_arr_if_not_len_enough(arr, len_required)
    return arr


def get_3_neighbouring_coord_values(coord_val, coord_resolution):
    """
    Get the three neighbouring values to a given input coord (coord_val)

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687

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
    if not isinstance(coord_val, float) or not isinstance(coord_resolution, float):
        coord_val = float(coord_val)
        coord_resolution = float(coord_resolution)

    return np.array(
        [coord_val - coord_resolution, coord_val, coord_val + coord_resolution]
    )


def quadratic_func(x, y):
    """
    Quadratic function.

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687

    Parameters
    ----------
    x : xr.DataArray or array-like
        Array 1
    y : xr.DataArray or array-like
        Array 2
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=np.RankWarning)
        p = np.polyfit(x, y, deg=2)
    return p


def apply_quadratic_func(x, y, vals):
    """
    Apply quadratic function to an array of values (vals).

    Component of method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    & Grise & Polvani (2016) https://doi.org/10.1002/2015JD024687


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

    data_ua = abs(data_row["ua"].dropna(dim="lat"))
    try:
        # Try and fit the parabola
        _, coeff = fit_parabola(data_ua["lat"].data, data_ua.data)
    except (RuntimeError, ValueError):
        coeff = [np.nan, np.nan]
    # Check that the max lat is not above the min or max (problem with using a parabola)
    if coeff[0] < data_row["lat"].min() or coeff[0] > data_row["lat"].max():
        coeff = [np.nan, np.nan]
    data_row["jet_lat"] = float(coeff[0])
    data_row["jet_speed"] = float(coeff[1])
    return data_row


def fit_parabola(x, y):
    """
    Fits a parabola.

    Component of method from Barnes & Polvani (2015) http://journals.ametsoc.org/doi/10.1175/JCLI-D-14-00589.1
    """
    coeff, _ = scipy.optimize.curve_fit(parabola, x, y)
    fitted_parabola = parabola(x, coeff[0], coeff[1], coeff[2])
    return fitted_parabola, coeff


def parabola(x, lat, speed, a):
    """
    Parabola.

    Component of method from Barnes & Polvani (2015) http://journals.ametsoc.org/doi/10.1175/JCLI-D-14-00589.1
    """
    return -a * (x - lat) ** 2 + speed


def assign_jet_lat_speed_to_ten_day_mean_data(ten_day_mean):
    """
    Joins the values for latitude and windspeed of the point of maximum windspeed for a given sector to the ten_day_mean data.

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
    Runs a cubic spline interpolation using 2 equal-length arrays.

    Component of method from Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1

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
    return scipy.interpolate.interp1d(x, y, kind="cubic", fill_value="extrapolate")


def run_cubic_spline_interpolation_to_get_max_lat_and_ws(
    data, lat_resolution, ua_col="ua"
):
    """
    Runs a cubic spline interpolation to find maximum latitude and maximum windspeed at a given resolution of latitude.

    Component of method from Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1

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
    scaled_lats = data_utils.rescale_lat_resolution(data["lat"], lat_resolution)
    csi = cubic_spline_interpolation(data["lat"], data[ua_col])
    interpolated_ws = csi(scaled_lats)
    max_lat = scaled_lats[np.argmax(interpolated_ws)]
    max_ws = max(interpolated_ws)
    return max_lat, max_ws


def run_cubic_spline_interpolation_for_each_unit_of_climatology_to_get_max_lat_and_ws(
    data, lat_resolution, time_col
):
    """
    Runs a cubic spline interpolation to find maximum latitude and maximum windspeed at a given resolution of latitude for each.

    Component of method from Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1

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
    seasonal_climatology = jetstream_metrics_components.get_climatology(data, "season")
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
    Will get the centroid latitude of the u-component wind using:
    jet_lat = integral(60deg, 30deg)(zonal_u**2*lat) dlat / integral(60deg, 30deg)(zonal_u**2) dlat

    Component of method from Ceppi et al (2018) https://doi.org/10.1175/JCLI-D-17-0323.1 & Zappa et al. 2018 https://doi.org/10.1029/2019GL083653

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


def get_latitude_value_in_data_row(data_row, lat_val):
    """
    Will return the latitude value of a given data row. Needed to loop through data and external array

    Component of method from Ceppi et al (2018) https://doi.org/10.1175/JCLI-D-17-0323.1 & Zappa et al. 2018 https://doi.org/10.1029/2019GL083653

    Parameters
    ----------
    data_row : xarray.DataArray
        Data row to extract variable from using given latitude
    lat_val : float or int
        latitude value to extract from data row
    """
    return float(data_row.sel(lat=lat_val, method="nearest"))


def get_moving_averaged_smoothed_jet_lats_for_one_day(data_row, width_of_pulse=10):
    """
    Get moving average of smoothed_jet_lats for one day. Used in combination with a groupby mapping.

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
    data_row["jet_lat"] = get_jet_lat_by_lon(data_row["ua"])
    data_row["smoothed_jet_lat"] = smooth_jet_lat_across_lon_with_rectangular_pulse(
        data_row["jet_lat"], width_of_pulse=width_of_pulse
    )
    return data_row


def smooth_jet_lat_across_lon_with_rectangular_pulse(jet_lat_data, width_of_pulse):
    """
    Smooth jet position (jet latitude) by carrying out a convolution with a rectangular pulse.

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
    Gets all the latitudes of maximum wind-speeds by longitude.

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
    argmax_lat_indices = xr.DataArray(
        dims=("lon",),
        data=np.nanargmax(data_row, axis=data_row.get_axis_num("lat")),
    )
    return data_row["lat"][argmax_lat_indices]
