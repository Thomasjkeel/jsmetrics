# -*- coding: utf-8 -*-

"""
    Jet-stream algorithms used in the literature.
    Algorithms are treated seperately to metrics, as in general, metrics are used to summarise information about the jet-stream, the algorithms simply identify it

    Classes and Functions ordered by paper publish year.
"""

# imports
from . import jet_core_algorithms_components
from jsmetrics.utils import windspeed_utils
from jsmetrics.core.check_data import sort_xarray_data_coords

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


@sort_xarray_data_coords(coords=["lat", "lon"])
def koch_et_al_2006(data, ws_threshold=30):
    """
    Calculates the weighted average windspeed and applies a threshold to identify the jet.
    The actual methodology uses 100-400 hPa and 30 ms^-1 as the windspeed threshold.

    weighted average windspeed = 1/(p2-p1) integral(p2, p1)(u^2+v^2)^(1/2)dp
    where p1, p2 is min, max pressure level

    Method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind
    ws_threshold : int or float
        Windspeed threshold for jet-stream (default: 30 ms-1)

    Returns
    ----------
    weighted_average_ws : xarray.Dataset
        data containing weighted average ws above windspeed threshold

    """
    if data["plev"].count() < 2:
        raise ValueError(
            "Need at least two plevs to calculate weighted average windspeed"
        )

    # Step 1: get all pressure levels (hPa) as list
    all_plevs_hPa = jet_core_algorithms_components.get_all_hPa_list(data)

    # Step 2: get weighted sum windspeed
    sum_weighted_ws = jet_core_algorithms_components.get_sum_weighted_ws(
        data, all_plevs_hPa
    )

    # Step 3: calculate average weighted
    weighted_average_ws = jet_core_algorithms_components.get_weighted_average_ws(
        sum_weighted_ws, all_plevs_hPa
    )

    # Step 4: Apply windspeed threshold
    weighted_average_ws = weighted_average_ws.where(weighted_average_ws >= ws_threshold)

    weighted_average_ws = weighted_average_ws.fillna(0.0)
    # Step 5: turn into dataset
    weighted_average_ws = weighted_average_ws.rename("weighted_average_ws").to_dataset()
    return weighted_average_ws


@sort_xarray_data_coords(coords=["lat", "lon"])
def schiemann_et_al_2009(data):
    """
    An occurrence-based jet climatology. Uses three rules: wind-vector wind-speed is local maxima is above 30 m/s and u-wind is more than 0 m/s.

    Method from Schiemann et al 2009 https://doi.org/10.1175/2008JCLI2625.1

    NOTE: Currently takes a very long time i.e. 8 seconds per time unit (i.e. 8 seconds per day) on AMD Ryzen 5 3600 6-core processor
    TODO: speed this metric up

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind

    Returns
    ----------
    output : xr.Dataset
        Data with local jet maximas
    """
    #  Step 1. Calculate wind vector
    data["ws"] = windspeed_utils.get_resultant_wind(data["ua"], data["va"])

    #  Step 2. Calculate jet maximas
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        output = jet_core_algorithms_components.get_local_jet_maximas_by_oneday_by_plev(
            data
        )
    else:
        output = data.groupby("time").map(
            jet_core_algorithms_components.get_local_jet_maximas_by_oneday_by_plev
        )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def manney_et_al_2011(data, ws_core_threshold=40, ws_boundary_threshold=30):
    """
    Looks to get seperate jet cores based on boundary and threshold. Core are discovered where 8-cells are above boundary threshold
    Paper uses 100-400 hPa.

    Method from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
    Also see Manney et al. 2011, 2014, 2017 and Manney & Hegglin 2018

    NOTE: Currently takes a long time i.e. 7.6 seconds per time unit with 8 plevs (i.e. 7.6 seconds per day) on AMD Ryzen 5 3600 6-core processor

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind
    ws_core_threshold : int or float
        Threshold used for jet-stream core point (default=40)
    ws_boundary_threshold : int or float
        Threshold for jet-stream boundary point (default=30)

    Returns
    ----------
    output : xarray.Dataset
        Data containing jet-cores (ID number relates to each unique core)
    """
    if "plev" not in data.dims:
        data = data.expand_dims("plev")

    # Step 1. Run Jet-stream Core Idenfication Algorithm
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        if "time" in data.dims:
            data = data.isel(time=0)
        output = jet_core_algorithms_components.run_jet_core_algorithm_on_one_day(
            data, ws_core_threshold, ws_boundary_threshold
        )
    else:
        output = data.groupby("time").map(
            jet_core_algorithms_components.run_jet_core_algorithm_on_one_day,
            (
                ws_core_threshold,
                ws_boundary_threshold,
            ),
        )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def penaortiz_et_al_2013(data):
    """
    Will calculate local wind maxima days per monthyear

    Method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    NOTE: Currently takes a long time i.e. 1.3 seconds per time unit with 8 plevs (i.e. 1.3 seconds per day) on AMD Ryzen 5 3600 6-core processor

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind

    Returns
    ----------
    output : xarray.Dataset
        Data containing number of days per month with local wind maxima
    """
    #  Step 1. Calculate wind vector
    data["ws"] = windspeed_utils.get_resultant_wind(data["ua"], data["va"])

    #  Step 2. Make array of zeros for local wind maxima location algorithm
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1 and "time" not in data.dims:
        data = data.expand_dims("time")
    output = jet_core_algorithms_components.get_empty_local_wind_maxima_data(data)

    #  Step 3. Find local wind maxima locations by day
    output = output.groupby("time").map(
        jet_core_algorithms_components.get_local_wind_maxima_by_timeunit
    )

    #  Step 4. Get number of days per month with local wind maxima
    output = jet_core_algorithms_components.get_number_of_timeunits_per_monthyear_with_local_wind_maxima(
        output
    )

    #  Step 5. Sort monthyear data into PJ and STJ
    output = jet_core_algorithms_components.subdivide_local_wind_maxima_into_stj_pfj(
        output
    )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def kuang_et_al_2014(data, occurence_ws_threshold=30):
    """
    Looks to get event-based jet occurrence and jet center occurrence of JS (1 is occurence, 2 is core).

    Method from Kuang et al (2014) https://doi.org/10.1007/s00704-013-0994-x

    NOTE: Currently takes a long time i.e. 2 seconds per time unit with 1 plev (i.e. 2 seconds per day) on AMD Ryzen 5 3600 6-core processor

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind
    occurence_ws_threshold : int or float
        Threshold used to identify a jet-stream occurence point (default=30)

    Returns
    ----------
    output : xarray.Dataset
        Data containing jet-occurence and jet-centres (1 is occurence, 2 is core)
    """
    if "plev" in data.dims:
        if data["plev"].count() == 1:
            data = data.isel(plev=0)
        else:
            raise ValueError("Please subset to one plev value for algorithm")

    # Step 1. Run Jet-stream Occurence and Centre Algorithm
    if "time" not in data.coords:
        output = (
            jet_core_algorithms_components.run_jet_occurence_and_centre_alg_on_one_day(
                data, occurence_ws_threshold
            )
        )
    else:
        if data["time"].size == 1:
            output = jet_core_algorithms_components.run_jet_occurence_and_centre_alg_on_one_day(
                data, occurence_ws_threshold
            )
        else:
            output = data.groupby("time").map(
                jet_core_algorithms_components.run_jet_occurence_and_centre_alg_on_one_day,
                (occurence_ws_threshold,),
            )
    return output
