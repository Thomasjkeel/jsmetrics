# -*- coding: utf-8 -*-

"""
    Metrics (or Indices) used to identify or classify jet-stream in the literature. This includes metrics to calculate
    latitudes at which highest wind-speed occurs or and upper-level wind sinuosity with the specific purpose of capturing the jet-stream 'waviness'

    All functions should return a xarray.Dataset.

    Classes and Functions ordered by paper publish year.
"""

# imports
import numpy as np
import xarray
from . import jet_statistics_components
from jsmetrics.utils import data_utils, spatial_utils, windspeed_utils
from jsmetrics.core.check_data import sort_xarray_data_coords

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


@sort_xarray_data_coords(coords=["lat", "lon"])
def archer_caldeira_2008(data):
    """
    Calculates the mass-weighted average wind speed, mass flux weighted pressure
    and mass flux weighted latitude. This method has some similarities to method
    used in Koch et al. 2006. In paper, 100-400 hPa is used.

    Method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind

    Returns
    ----------
    output : xarray.Dataset
        Data containing mass weighted average ws, mass flux weighted pressure and latitude
    """
    #  Step 1. Get monthly means
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        print(
            "Warning: only found one time step, and the time coord is being renamed month."
        )
        mon_mean = data.assign_coords(month=data["time"].dt.month)
        mon_mean = mon_mean.expand_dims("month")
    else:
        mon_mean = data.groupby("time.month").mean()

    #  Step 2. Calculate wind-speed from u and v-component wind
    mon_mean["ws"] = windspeed_utils.get_resultant_wind(mon_mean["ua"], mon_mean["va"])

    #  Step 3. Calculate mass weighted average
    mass_weighted_average = jet_statistics_components.calc_mass_weighted_average(
        mon_mean, ws_col="ws"
    )
    mass_flux_weighted_pressure = (
        jet_statistics_components.calc_mass_flux_weighted_pressure(
            mon_mean, ws_col="ws"
        )
    )
    mass_flux_weighted_latitude = (
        jet_statistics_components.calc_mass_flux_weighted_latitude(
            mon_mean, lat_min=15, lat_max=75, ws_col="ws"
        )
    )
    output = data.assign(
        {
            "mass_weighted_average_ws": (
                ("month", "lat", "lon"),
                mass_weighted_average.data,
            ),
            "mass_flux_weighted_pressure": (
                ("month", "lat", "lon"),
                mass_flux_weighted_pressure.data,
            ),
            "mass_flux_weighted_latitude": (
                ("month", "lon"),
                mass_flux_weighted_latitude.data,
            ),
        }
    )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def woollings_et_al_2010(data, filter_freq=10, window_size=61):
    """
    Follows an in-text description of 4-steps describing the algorithm of jet-stream identification from Woollings et al. (2010).
    Will calculate this metric based on data (regardless of pressure level of time span etc.).

    Method from Woollings et al (2010) http://dx.doi.org/10.1002/qj.625

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- component wind
    filter_freq : int
        number of days in filter (default=10 days)
    window_size : int
        number of days in window for Lancoz filter (default=61 days)

    Returns
    ----------
    fourier_filtered_data : xarray.Dataset
        Data containing maximum latitudes and maximum windspeed at those lats and fourier-filtered versions of those two variables
    """
    if isinstance(data, xarray.DataArray):
        data = data.to_dataset()
    # Step 1: Calculate long and/or plev mean
    zonal_mean = windspeed_utils.get_zonal_mean(data)

    # Step 2: Apply n-day lancoz filter
    lancoz_filtered_mean_data = jet_statistics_components.apply_lanczos_filter(
        zonal_mean["ua"], filter_freq, window_size
    )

    # Step 3: Calculate max windspeed and lat where max ws found
    max_lat_ws = np.array(
        list(
            map(
                jet_statistics_components.get_latitude_and_speed_where_max_ws,
                lancoz_filtered_mean_data[:],
            )
        )
    )
    zonal_mean_lat_ws = jet_statistics_components.assign_jet_lat_and_speed_to_data(
        zonal_mean, max_lat_ws
    )
    # Step 4: Make seasonal climatology
    climatology = jet_statistics_components.get_climatology(zonal_mean_lat_ws, "season")

    # Step 5: Apply low-freq fourier filter to both max lats and max ws
    fourier_filtered_lats = jet_statistics_components.apply_low_freq_fourier_filter(
        climatology["jet_lat"].values, highest_freq_to_keep=2
    )
    fourier_filtered_ws = jet_statistics_components.apply_low_freq_fourier_filter(
        climatology["jet_speed"].values, highest_freq_to_keep=2
    )
    time_dim = climatology["jet_speed"].dims[0]
    output = jet_statistics_components.assign_filtered_lats_and_ws_to_data(
        zonal_mean_lat_ws,
        fourier_filtered_lats.real,
        fourier_filtered_ws.real,
        dim=time_dim,
    )

    # Step 6: Calculate jet latitude and jet speed anomalies from the seasonal cycle

    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def barnes_polvani_2013(data, filter_freq=10, window_size=41):
    """
    Pressure weighted u-component wind then gets low-pass lanczos filtered (10-day, 41 weights) and 0.01 quadratic function applied
    for jet-lat and speed. "We define the jet width as the full width at half of the maximum jet speed".

    Method from Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component wind
    filter_freq : int
        number of days in filter (default=10 timeunits)
    window_size : int
        number of days in window for low-pass Lancoz filter (default=41 timeunits)

    Returns
    ----------
    output : xarray.Dataset
        Data containing values for z_lat, z_spd, z_width for jet-stream latitude, speed and width
    """
    #  Step 1. Get pressure-weighted u-component wind
    pressure_weighted_ua = jet_statistics_components.calc_mass_weighted_average(
        data, ws_col="ua"
    )
    #  Step 2. Filter pressure-weighted u-component wind with low-pass Lanczos filter
    filtered_pressure_weighted_ua = jet_statistics_components.apply_lanczos_filter(
        pressure_weighted_ua,
        filter_freq=filter_freq,
        window_size=window_size,
    )
    #  Step 3. Turn dataarray into dataset for next part
    filtered_pressure_weighted_ua = filtered_pressure_weighted_ua.to_dataset(name="ua")

    #  Step 4.  Get max latitude and wind speed at max
    zonal_mean = windspeed_utils.get_zonal_mean(filtered_pressure_weighted_ua)
    all_max_lats_and_ws = np.array(
        list(
            map(
                jet_statistics_components.get_3_latitudes_and_speed_around_max_ws,
                zonal_mean["ua"],
            )
        )
    )

    #  Step 5. Scale max latitude to 0.01 degree using quadratic function
    scaled_max_lats = []
    scaled_max_ws = []
    for max_lat_and_ws in all_max_lats_and_ws:
        if not np.isnan(np.min(max_lat_and_ws)):
            (
                scaled_max_lat,
                scaled_ws,
            ) = jet_statistics_components.get_latitude_and_speed_where_max_ws_at_reduced_resolution(
                max_lat_and_ws, lat_resolution=0.01
            )
            scaled_max_lats.append(scaled_max_lat)
            scaled_max_ws.append(scaled_ws)
        else:
            if np.isnan(max_lat_and_ws).all():
                scaled_max_lats.append(np.nan)
                scaled_max_ws.append(np.nan)
            else:
                if np.nanmax(max_lat_and_ws[0]) == data["lat"].max():
                    scaled_max_lats.append(np.nanmax(max_lat_and_ws[0]))
                    scaled_max_ws.append(np.nanmax(max_lat_and_ws[1]))
                elif np.nanmin(max_lat_and_ws[0]) == data["lat"].min():
                    scaled_max_lats.append(np.nanmin(max_lat_and_ws[0]))
                    scaled_max_ws.append(np.nanmin(max_lat_and_ws[1]))
                else:
                    scaled_max_lats.append(np.nan)
                    scaled_max_ws.append(np.nan)

    #  Step 6. Get jet-widths using scaled windspeed and usual jet-lat
    max_lats = all_max_lats_and_ws[::, 0, 1]
    jet_widths = list(
        map(
            lambda zm, la, wa: jet_statistics_components.calc_jet_width_for_one_day(
                zm, la, wa
            ),
            zonal_mean["ua"],
            max_lats,
            scaled_max_ws,
        )
    )

    output = data.assign(
        {
            "jet_lat": (("time"), scaled_max_lats),
            "jet_speed": (("time"), scaled_max_ws),
            "jet_width": (("time"), jet_widths),
        }
    )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def barnes_polvani_2015(data):
    """
    Calculates the jet speed and jet position by fitting a parabola around the
    maximum of zonally average wind and taking the maximum magnitude and position
    to be the jet speed and jet latitude respectively.

    Method from Barnes & Polvani (2015) https://doi.org/10.1175/JCLI-D-14-00589.1

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component wind

    Returns
    ----------
    output : xarray.Dataset
        Data containing jet-stream position and jet-speed
    """
    # Step 1. Get zonal mean
    zonal_mean = windspeed_utils.get_zonal_mean(data)

    # Step 2. Get jet lat and jet speed values
    if zonal_mean["time"].size == 1:
        output = jet_statistics_components.get_jet_lat_and_speed_using_parabola_by_day(
            zonal_mean
        )
    else:
        output = zonal_mean.groupby("time").map(
            jet_statistics_components.get_jet_lat_and_speed_using_parabola_by_day
        )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def barnes_simpson_2017(data):
    """
    "Time series of jet latitude and jet speed are defined as the latitude and speed of the 10-day-averaged
     maximum 700-hPa zonal winds averaged over the longitudinal sector of interest"

     Method from Barnes & Simpson 2017 https://doi.org/10.1175/JCLI-D-17-0299.1

     Parameters
     ----------
     data : xarray.Dataset
         Data containing u-component wind

     Returns
     ----------
     output : xarray.Dataset
         Data with max latitude and max windspeed for North Atlantic (280.E to 350. E) and North Pacific (120.E to 230. E) sectors
    """
    if "plev" in data.dims:
        if data["plev"].count() == 1:
            data = data.isel(plev=0)
        else:
            print(
                "this metric was meant to only work on one plev, please subset plev to one value. For now taking the mean..."
            )
            data = data.mean("plev")
    data = data.mean("lon")
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1 and "time" not in data.dims:
        data = data.expand_dims("time")
    if not data.indexes["time"].is_monotonic_increasing:
        raise IndexError("Data needs to have a montonic increasing index")
    # Check that data can be resampled into 10 days
    if not data["time"].size == 1:
        time_step_in_data = int((data["time"][1] - data["time"][0]).dt.days)
        if time_step_in_data <= 10:
            data = data.resample(time="10D").mean()
            time_step_in_data = 10
        else:
            print(
                f"Warning this method was developed for 10 day average and data has larger time-step than 10 days. Time step is {time_step_in_data} days"
            )
    #  Drop all NaN slices
    data = data.dropna("time")
    data = jet_statistics_components.calc_latitude_and_speed_where_max_ws(data)
    return data


@sort_xarray_data_coords(coords=["lat", "lon"])
def grise_polvani_2017(data):
    """
    Calculates maximum latitude of jet-stream to 0.01 degree resolution each time unit
    Method from Grise & Polvani (2017) https://doi.org/10.1175/JCLI-D-16-0849.1

    See also Ceppi et al. 2012
    Methodology is for Southern Hemisphere
    NOTE: This method also uses poleward edge of sub-tropical dry zone and poleward edge of Hadley cell derived from precip. record

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component wind

    Returns
    ----------
    output : xarray.Dataset
        Data containing max latitudes per time unit scaled to 0.01 resolution
    """
    if isinstance(data, xarray.DataArray):
        data = data.to_dataset()

    # Step 0: Expand time dimensions so we can map a function to the dataset properly
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        data = data.expand_dims("time")

    # Step 1. Calculate zonal-mean
    zonal_mean = windspeed_utils.get_zonal_mean(data)

    # Step 2. Get the 3 latitudes and speeds around max zonal wind-speed (e.g. lat-1, lat, lat+1)
    all_max_lats_and_ws = np.array(
        list(
            map(
                jet_statistics_components.get_3_latitudes_and_speed_around_max_ws,
                zonal_mean["ua"],
            )
        )
    )

    #  Step 3. Apply quadratic function to get max latitude at 0.01 degree resolution
    scaled_max_lats = []
    scaled_max_ws = []
    for max_lat_and_ws in all_max_lats_and_ws:
        try:
            (
                scaled_max_lat,
                scaled_ws,
            ) = jet_statistics_components.get_latitude_and_speed_where_max_ws_at_reduced_resolution(
                max_lat_and_ws, lat_resolution=0.01
            )
        except Exception as e:
            print(e)
            scaled_max_lat = np.nan
            scaled_ws = np.nan
        scaled_max_lats.append(scaled_max_lat)
        scaled_max_ws.append(scaled_ws)

    #  Step 4. Assign scaled max lats back to data
    output = data.assign(
        {
            "jet_lat": (("time"), scaled_max_lats),
            "jet_speed": (("time"), scaled_max_ws),
        }
    )
    return output


# @sort_xarray_data_coords(coords=["lat", "lon"])
# def molnos_et_al_2017(data):
#     """
#     Write function description
#     """
#     return data


@sort_xarray_data_coords(coords=["lat", "lon"])
def bracegirdle_et_al_2018(data):
    """
    Calculates the seasonal and annual jet-stream position from a cubic spline interpolation of zonal wind climatology.
    Method from Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1

    NOTE: Originally for Southern Hemisphere

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component windspeed

    Returns
    ----------
    output : xarray.Dataset
        Data containing seasonal and annual jet-stream position and strength (ms-1)
    """
    if isinstance(data, xarray.DataArray):
        data = data.to_dataset()
    if "plev" in data.dims:
        if data["plev"].count() == 1:
            data = data.isel(plev=0)
        else:
            print(
                "this metric was meant to only work on one plev, please subset plev to one value. For now taking the mean..."
            )
            data = data.mean("plev")
            # raise ValueError("Please subset to one plev value for this metric")

    # Step 0: Expand time dimensions so we can map a function to the dataset properly
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        data = data.expand_dims("time")

    #  Step 1. Make seasonal & annual climatologies
    seasonal_climatology = jet_statistics_components.get_climatology(data, "season")
    annual_climatology = jet_statistics_components.get_climatology(data, "year")

    #  Step 2. Get zonal mean from climatologies
    seasonal_zonal_mean = seasonal_climatology.mean("lon")
    annual_zonal_mean = annual_climatology.mean("lon")

    #  Step 3. Cubic spline interpolation to each climatology at latitude resolution of 0.075 degrees
    (
        seasonal_max_lats,
        seasonal_max_ws,
    ) = jet_statistics_components.run_cubic_spline_interpolation_for_each_unit_of_climatology_to_get_max_lat_and_ws(
        seasonal_zonal_mean, lat_resolution=0.075, time_col="season"
    )
    (
        annual_max_lats,
        annual_max_ws,
    ) = jet_statistics_components.run_cubic_spline_interpolation_for_each_unit_of_climatology_to_get_max_lat_and_ws(
        annual_zonal_mean, lat_resolution=0.075, time_col="year"
    )

    # Step 4. Assign jet-stream position (JPOS) and jet-stream strength (JSTR) back to data
    output = data.assign(
        {
            "seasonal_JPOS": (("season"), seasonal_max_lats),
            "annual_JPOS": (("year"), annual_max_lats),
            "seasonal_JSTR": (("season"), seasonal_max_ws),
            "annual_JSTR": (("year"), annual_max_ws),
        }
    )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def ceppi_et_al_2018(data, lon_resolution=None):
    """
    Calculates the jet latitude per time unit where jet-lat is defined as a centroid of a zonal wind distribution.
    This method has been slightly adapted to include a jet speed extraction (after Screen et al. 2022 and refs therein).
    Method from Ceppi et al (2018) https://doi.org/10.1175/JCLI-D-17-0323.1

    Also see Zappa et al. 2018 method which includes exclusion of <0 m/s u-wind

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component windspeed
    lon_resolution : numeric
        Resolution to use for longitude coord if size 1

    Returns
    ----------
    output : xarray.Dataset
        Data containing centroid latitude of u-wind for each time unit (e.g. each day)
    """
    #  Step 1. Get area in m2 by latitude/longitude grid cells
    if not data["lon"].size == 1 and not data["lat"].size == 1:
        total_area_m2 = spatial_utils.grid_cell_areas(data["lon"], data["lat"])
    elif lon_resolution and not data["lat"].size == 1 and data["lon"].size == 1:
        lons_to_use = [float(data["lon"]), float(data["lon"]) + lon_resolution]
        total_area_m2 = spatial_utils.grid_cell_areas(lons_to_use, data["lat"])
        total_area_m2 = total_area_m2.mean(axis=1)
        total_area_m2 = total_area_m2.reshape(-1, 1)
        if data["lon"].shape == ():
            data = data.expand_dims("lon")
    else:
        raise ValueError(
            "For this method, your data needs to have at least 2 'lat' values and 'lon' values needs to be more than one unless you set the 'lon_resolution' parameter"
        )

    data["total_area_m2"] = (("lat", "lon"), total_area_m2)

    #  Step 2. calculate zonal mean
    zonal_mean = windspeed_utils.get_zonal_mean(data)

    # Step 3: Assign laitude of jet-stream centroids to main data
    data["jet_lat"] = jet_statistics_components.calc_centroid_jet_lat_from_zonal_mean(
        zonal_mean, area_by_lat=zonal_mean["total_area_m2"]
    )

    # Expand time dimension
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        data = data.expand_dims("time")
        zonal_mean = zonal_mean.expand_dims("time")

    # Step 4 (adapted from original methodology): Get nearest latitude actually in data to the one determined by metric
    nearest_latitudes_to_jet_lat_estimates = np.array(
        list(
            map(
                lambda row: data_utils.find_nearest_value_to_array(
                    data["lat"], float(row)
                ),
                data["jet_lat"],
            )
        )
    )

    # Step 5 (adapted from original methodology): Get speed of associated nearest latitude
    data["jet_speed"] = (
        ("time",),
        np.array(
            list(
                map(
                    lambda data_row, lat_val: jet_statistics_components.get_latitude_value_in_data_row(
                        data_row, lat_val
                    ),
                    zonal_mean["ua"],
                    nearest_latitudes_to_jet_lat_estimates,
                )
            )
        ),
    )
    return data


def zappa_et_al_2018(data, lon_resolution=None):
    """
    Calculates the jet latitude per time unit where jet-lat is defined as a centroid of a zonal wind distribution.
    This method has been slightly adapted to include a jet speed extraction (after Screen et al. 2022 and refs therein).
    Method from Zappa et al. 2018 https://doi.org/10.1029/2019GL083653
    Adapted from and very similar to Ceppi et al (2018).

    Also used in Ayres & Screen, 2019 and Screen et al. 2022. Similar methods used in: Chen et al. 2008; Ceppi et al. 2014, Ceppi et al. 2018.

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component windspeed
    lon_resolution : numeric
        Resolution to use for longitude coord if size 1

    Returns
    ----------
    output : xarray.Dataset
        Data containing centroid latitude of u-wind for each time unit (e.g. each day)
    """
    #  Step 1. Get area in m2 by latitude/longitude grid cells
    if not data["lon"].size == 1 and not data["lat"].size == 1:
        total_area_m2 = spatial_utils.grid_cell_areas(data["lon"], data["lat"])
    elif lon_resolution and not data["lat"].size == 1 and data["lon"].size == 1:
        lons_to_use = [float(data["lon"]), float(data["lon"]) + lon_resolution]
        total_area_m2 = spatial_utils.grid_cell_areas(lons_to_use, data["lat"])
        total_area_m2 = total_area_m2.mean(axis=1)
        total_area_m2 = total_area_m2.reshape(-1, 1)
        if data["lon"].shape == ():
            data = data.expand_dims("lon")
    else:
        raise ValueError(
            "For this method, your data needs to have at least 2 'lat' values and 'lon' values needs to be more than one unless you set the 'lon_resolution' parameter"
        )

    data["total_area_m2"] = (("lat", "lon"), total_area_m2)

    #  Step 2. calculate zonal mean and floor ua values to 0
    zonal_mean = windspeed_utils.get_zonal_mean(data)
    zonal_mean["ua"] = zonal_mean["ua"].where(lambda x: x > 0, 0)

    # Step 3: Assign laitude of jet-stream centroids to main data
    data["jet_lat"] = jet_statistics_components.calc_centroid_jet_lat_from_zonal_mean(
        zonal_mean, area_by_lat=zonal_mean["total_area_m2"]
    )

    # Expand time dimension
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        data = data.expand_dims("time")
        zonal_mean = zonal_mean.expand_dims("time")

    # Step 4 (adapted from original methodology): Get nearest latitude actually in data to the one determined by metric
    nearest_latitudes_to_jet_lat_estimates = np.array(
        list(
            map(
                lambda row: data_utils.find_nearest_value_to_array(
                    data["lat"], float(row)
                ),
                data["jet_lat"],
            )
        )
    )

    # Step 5 (adapted from original methodology): Get speed of associated nearest latitude
    data["jet_speed"] = (
        ("time",),
        np.array(
            list(
                map(
                    lambda data_row, lat_val: jet_statistics_components.get_latitude_value_in_data_row(
                        data_row, lat_val
                    ),
                    zonal_mean["ua"],
                    nearest_latitudes_to_jet_lat_estimates,
                )
            )
        ),
    )
    return data


@sort_xarray_data_coords(coords=["lat", "lon"])
def kerr_et_al_2020(data, width_of_pulse=10):
    """
    Described in section 2.4.2 of paper. Defines the latitude of the jet-stream as where the
    maximum zonal winds occur for each longitude for each time unit (i.e. day) before smoothing
    with a rectangular pulse (of width 10 degrees) to get a moving average.
    Method from Kerr et al. (2020) https://onlinelibrary.wiley.com/doi/10.1029/2020JD032735

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component windspeed at one plev

    Returns
    ----------
    output : xarray.Dataset
        Data containing jet-stream latitude by longitude and smoothed jet_latitude
    """
    if "plev" in data.dims:
        if data["plev"].count() == 1:
            data = data.isel(plev=0)
        else:
            print(
                "this metric was meant to only work on one plev, please subset plev to one value. For now taking the mean..."
            )
            data = data.mean("plev")
            # raise ValueError("Please subset to one plev value for this metric")
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    elif data["time"].size == 1:
        output = (
            jet_statistics_components.get_moving_averaged_smoothed_jet_lats_for_one_day(
                data, width_of_pulse
            )
        )
    else:
        output = data.groupby("time").map(
            jet_statistics_components.get_moving_averaged_smoothed_jet_lats_for_one_day,
            (width_of_pulse,),
        )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def blackport_fyfe_2022(data):
    """
    Described in MATERIALS AND METHODS section of this paper under: 'Metrics and analysis':
    The speed and latitude of the North Atlantic jet stream are calculated similar to a previous analysis (Barnes & Polvani 2015),
    with some minor differences. At each longitude between 60°W and 0°, we fit a parabola around the seasonal averaged,
    maximum U700 speeds found between 15°N and 75°N. The parabola is fit around the two points on either side of the maximum (five points in total).
    The maximum of the resulting function is the jet speed, and the latitude of this maximum is the jet latitude. If the maximum speed occurs
    at the northern or southern edge, then the speed and latitude are ignored for that longitude.

    Method from Blackport & Fyfe (2022) https://www.science.org/doi/10.1126/sciadv.abn3112

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component windspeed at one plev

    Returns
    ----------
    output : xarray.Dataset
        Data containing jet-stream latitude by longitude and smoothed jet_latitude
    """
    output = data
    return output
