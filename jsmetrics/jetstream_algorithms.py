# -*- coding: utf-8 -*-

"""
    Jet-stream algorithms used in the literature.
    Algorithms are different from metrics because ...
"""

# imports
# import numpy as np
# import xarray
from . import jetstream_metrics_utils
from . import general_utils, windspeed_utils

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def koch_et_al_2006(data, ws_threshold=30):
    """
    Method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    Calculates the weighted average windspeed and applies a threshold to identify the jet.
    The actual methodology uses 100-400 hPa and 30 ms^-1 as the windspeed threshold.

    weighted average windspeed = 1/(p2-p1) integral(p2, p1)(u^2+v^2)^(1/2)dp
    where p1, p2 is min, max pressure level

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

    TODO: check with chris
    TODO: add equation to this doc
    TODO: what if mbar?
    """
    if data["plev"].count() < 2:
        raise ValueError(
            "Need at least two plevs to calculate weighted average windspeed"
        )

    # Step 1: get all pressure levels (hPa) as list
    all_plevs_hPa = general_utils.get_all_hPa_list(data)

    # Step 2: get weighted sum windspeed
    sum_weighted_ws = jetstream_metrics_utils.get_sum_weighted_ws(
        data, all_plevs_hPa
    )

    # Step 3: calculate average weighted
    weighted_average_ws = jetstream_metrics_utils.get_weighted_average_ws(
        sum_weighted_ws, all_plevs_hPa
    )

    # Step 4: Apply windspeed threshold
    weighted_average_ws = weighted_average_ws.where(
        weighted_average_ws >= ws_threshold
    )

    weighted_average_ws = weighted_average_ws.fillna(0.0)
    # Step 5: turn into dataset
    weighted_average_ws = weighted_average_ws.rename(
        "weighted_average_ws"
    ).to_dataset()
    return weighted_average_ws


def archer_caldeira_2008(data):
    """
    Method from Archer & Caldiera (2008) https://doi.org/10.1029/2008GL033614

    Calculates the mass-weighted average wind speed, mass flux weighted pressure
    and mass flux weighted latitude. This method has some similarities to method
    used in Koch et al. 2006. In paper, 100-400 hPa is used.

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
    mon_mean = data.groupby("time.month").mean()

    #  Step 2. Calculate wind-speed from u and v-component wind
    mon_mean["ws"] = windspeed_utils.get_resultant_wind(
        mon_mean["ua"], mon_mean["va"]
    )

    #  Step 3. Calculate mass weighted average
    mass_weighted_average = jetstream_metrics_utils.calc_mass_weighted_average(
        mon_mean, ws_col="ws"
    )
    mass_flux_weighted_pressure = (
        jetstream_metrics_utils.calc_mass_flux_weighted_pressure(
            mon_mean, ws_col="ws"
        )
    )
    mass_flux_weighted_latitude = (
        jetstream_metrics_utils.calc_mass_flux_weighted_latitude(
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


def schiemann_et_al_2009(data):
    """
    Method from Schiemann et al 2009 https://doi.org/10.1175/2008JCLI2625.1

    Actual methodology uses 100-500 hPa
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
    output = data.groupby("time").map(
        jetstream_metrics_utils.get_local_jet_maximas_by_timeunit_by_plev
    )
    return output


def manney_et_al_2011(data, ws_core_threshold=40, ws_boundary_threshold=30):
    """
    Method from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
    Also see Manney et al. 2011, 2014, 2017 and Manney & Hegglin 2018

    Looks to get seperate jet cores based on boundary and threshold. Core are discovered where 8-cells are above boundary threshold
    Paper uses 100-400 hPa.
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
    output = data.groupby("time").map(
        jetstream_metrics_utils.run_jet_core_algorithm_on_one_day,
        (
            ws_core_threshold,
            ws_boundary_threshold,
        ),
    )
    return output


def penaortiz_et_al_2013(data):
    """
    Method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Will calculate local wind maxima days per monthyear
    Actual methodology uses 100-400 hPa

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
    local_wind_maxima = (
        jetstream_metrics_utils.get_empty_local_wind_maxima_data(data)
    )

    #  Step 3. Find local wind maxima locations by day
    local_wind_maxima_by_timeunit = local_wind_maxima.groupby("time").map(
        jetstream_metrics_utils.get_local_wind_maxima_by_timeunit
    )

    #  Step 4. Get number of days per month with local wind maxima
    local_wind_maxima_timeunits_by_monthyear = jetstream_metrics_utils.get_number_of_timeunits_per_monthyear_with_local_wind_maxima(
        local_wind_maxima_by_timeunit
    )
    local_wind_maxima_timeunits_by_monthyear = (
        local_wind_maxima_timeunits_by_monthyear.to_dataset()
    )

    #  Step 5. Sort into PJ and STJ
    output = jetstream_metrics_utils.subdivide_local_wind_maxima_into_stj_pfj(
        local_wind_maxima_timeunits_by_monthyear
    )
    return output


def kuang_et_al_2014(data, occurence_ws_threshold=30):
    """
    Method from Kuang et al (2014) https://doi.org/10.1007/s00704-013-0994-x

    Looks to get event-based jet occurrence and jet center occurrence of JS (1 is occurence, 2 is core).
    Best at 100-500 hPa
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
    output = data.groupby("time").map(
        jetstream_metrics_utils.run_jet_occurence_and_centre_alg_on_one_day,
        (occurence_ws_threshold,),
    )
    return output


# def kern_et_al_2018(data):
#     """
#     Write function description
#     TODO: ask about equation
#     """
#     return data
