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
    r"""
    This method follows a two-step procedure used to detect 'jet-event occurences'.
    The first step is to calculate the weighted average windspeed and then the
    second step is to apply a windspeed threshold to isolate jet events from that weighted average.
    The original methodology uses windspeed between 100-400 hPa to calculated the weighted average
    and 30 meters per second as the windspeed threshold.

    The weighted average windspeed for the jet events is calculated as follows:

    .. math::
        \alpha vel =  \frac{1}{p2-p1} \int_{p1}^{p2} (u^2+v^2)^{1/2} \,dp

    where p1, p2 is min, max pressure level.

    This method was first introduced in Koch et al (2006) (https://doi.org/10.1002/joc.1255)
    and is described in section 2.2.2 of that study. The original methodology provides a third step
    (to produce a climatology of jet events), but this has been ignored in this implementation.
    Instead, we have provided an example of how to calculate this after running this method
    in 'Examples' below.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    ws_threshold : int or float
        Windspeed threshold used to extract from weighted average (default: 30 ms-1)

    Returns
    ----------
    xarray.Dataset
        A dataset containing weighted average ws above windspeed threshold

    Notes
    -----
    This equation for this method is provided on pg 287 of the Koch et al. 2006 paper.
    In the original paper, they accumulate the jet events into two-class jet typology (described in section 2.2.3
    of Koch et al. 2006)

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        uv_data = xr.open_dataset('path_to_uv_data')

        # Subset dataset to range used in original methodology (100-400 hPa)):
        uv_sub = uv_data.sel(plev=slice(100, 400))

        # Run algorithm:
        koch_outputs = jsmetrics.jet_core_algorithms.koch_et_al_2006(uv_sub, ws_threshold=30)

        # Produce climatology of jet occurence events for each season and each month:
        koch_month_climatology = koch_outputs.groupby("time.month").mean("time")
        koch_season_climatology = koch_outputs.groupby("time.season").mean("time")

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

    # Step 4: Apply windspeed threshold to get jet event dataset
    jet_events = weighted_average_ws.where(weighted_average_ws >= ws_threshold)

    jet_events = jet_events.fillna(0.0)

    # Step 5: turn into dataset
    jet_event_ds = jet_events.rename("jet_events_ws").to_dataset()
    return jet_event_ds


@sort_xarray_data_coords(coords=["lat", "lon"])
def schiemann_et_al_2009(data, ws_threshold=30):
    r"""
    This method detects jet occurrences, whereby each jet occurence is detected based
    on three rules applied to inputted wind speed (V = [u, v]):
        1. \|V\| is a local maxima in latitude and altitude plane
        2. \|V\| ≥ 30 m/s
        3. \|u\| ≥ 0 m/s.

    This method was originally introduce in Schiemann et al 2009 (https://doi.org/10.1175/2008JCLI2625.1)
    and is described in Section 2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.


    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    ws_threshold : int or float
        Windspeed threshold used to extract jet maxima (default: 30 ms-1)

    Returns
    ----------
    output : xr.Dataset
        Data with local jet maximas

    Notes
    -----
    While the original method is built on a four dimension slice of wind speed (time, lat, lon, plev),
    This implementation will work where there is only one pressure level, so a 3-d slice (time, lat, lon).

    **Slow method:** Due to the nature of this method, it currently takes a very long time to run,
    i.e. 8 seconds per time unit on AMD Ryzen 5 3600 6-core processor.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        uv_data = xr.open_dataset('path_to_uv_data')

        # Subset dataset to range used in original methodology (100-500 hPa & 16.7-58.25 N, 42.5-220.5 E)):
        uv_sub = uv_data.sel(plev=slice(100, 500), lat=slice(16.7, 58.25), lon=slice(42.5, 220.5))

        # Run algorithm:
        schiemann_outputs = jsmetrics.jet_core_algorithms.schiemann_et_al_2009(uv_sub, ws_threshold=30)

        # Produce a jet occurence count across all pressure levels
        schiemann_jet_counts_all_levels = schiemann_outputs['jet_occurence'].sum(('time', 'plev'))

        # Use the jet occurence values as a mask to extract the jet windspeeds
        schiemann_jet_ws = schiemann_outputs.where(schiemann_outputs['jet_occurence'] > 0)['ws']
    """
    #  Step 1. Calculate wind vector
    data["ws"] = windspeed_utils.get_resultant_wind(data["ua"], data["va"])

    #  Step 2. Calculate jet occurences
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        if "time" in data.dims:
            data = data.squeeze("time")
        output = (
            jet_core_algorithms_components.get_local_jet_occurence_by_oneday_by_plev(
                data, ws_threshold=ws_threshold
            )
        )
    else:
        output = data.groupby("time").map(
            jet_core_algorithms_components.get_local_jet_occurence_by_oneday_by_plev,
            (ws_threshold,),
        )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def manney_et_al_2011(
    data,
    jet_core_plev_limit,
    jet_core_ws_threshold=40,
    jet_boundary_ws_threshold=30,
    ws_drop_threshold=25,
    jet_core_lat_distance=15,
):
    r"""
    This method detects jet cores and defines a boundary region beside those cores based on two windspeed thresholds.
    Two additional checks are applied after initial detection of cores to check whether cores within the same windspeed region
    are part of the same feature (default is 30 m/s, see 'jet_boundary_ws_threshold').
    These two checks are achieved by checking whether regions with multiple jet cores are more than a certain distance apart
    (default is 15 degrees, see 'jet_core_lat_distance') and the windspeed between two cores does not drop below a threshold
    (default is 25 m/s, see 'ws_drop_threshold')

    This method returns four outputs
        1. 'jet_core_mask' -- Regions within each latitude/altitude that are local maxima have windspeeds above the 'jet_core_ws_threshold'
        2. 'jet_region_mask' -- Regions above, below, left and right of the jet core with windspeed above the 'jet_boundary_ws_threshold'
        3. 'jet_region_above_ws_threshold_mask' -- All contigious regions of windspeeds emcompassing a jet core above the 'jet_boundary_ws_threshold' (i.e. not just above, below, left and right)
        4. 'ws' -- Wind speed calculated from 'ua', 'va' inputs.

    This method was originally introduce in Manney et al. (2011) (https://doi.org/10.5194/acp-11-6115-2011),
    and is described in Section 3.1 of that study. This method is also known as JETPAC, and available in its
    original form from NASA JPL.

    There is an update to this method introduced in Manney & Hegglin 2018 to include physically-based method to extract the
    subtropical jet is identified (and thus distinguished from polar jets).


    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    jet_core_plev_limit: tuple or array
        Sequence of two values relating to the pressure level limit of the jet cores (original paper uses 100hPa 400 hPa)
    jet_core_ws_threshold : int or float
        Threshold used for jet-stream core point (default=40 m/s)
    jet_boundary_ws_threshold : int or float
        Threshold for jet-stream boundary point (default=30 m/s)
    ws_drop_threshold : int or float
        Threshold for drop in windspeed along the line between cores (default: 25 m/s)
    jet_core_lat_distance : int or float
        Threshold for maximum distance between cores to be counted the same (default: 15 degrees)

    Returns
    ----------
    output : xarray.Dataset
        Data containing the variable 'jet-core_id' (ID number relates to each unique core)

    Notes
    -----
    The implementation of this method varies slightly from the original, in that this method will return
    variables that have 0, 1+ values, so that the user can use these as a mask on other variables such as windspeed
    (see 'Examples' for demonstration of how to use the mask).
    Also, 'jet_region_above_ws_threshold_mask' is provided here as a alternative to using a contour to check which regions
    encompass jet cores.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        uv_data = xr.open_dataset('path_to_uv_data')

        # Subset dataset to range appropriate for original methodology (100-1000 hPa)):
        uv_sub = uv_data.sel(plev=slice(100, 1900))

        # Run algorithm:
        manney_outputs = jsmetrics.jet_core_algorithms.manney_et_al_2011(uv_sub, ws_core_threshold=40, ws_boundary_threshold=30, jet_core_plev_limit=(100, 400))

        # Use the jet core mask to extract the jet windspeeds
        manney_jet_ws = manney_outputs.where(manney_outputs['jet_core_mask'])['ws']

    """
    if "plev" not in data.dims:
        data = data.expand_dims("plev")

    if not jet_core_plev_limit:
        raise KeyError(
            "Please provide a pressure level limit for jet cores returned by this metric. As an example the original methodology used 100-400 hPa as a limit (to replicate this, pass the parameter jet_core_plev_limit=(100, 400))"
        )

    # Step 1. Calculate wind speed from ua and va components.
    data["ws"] = windspeed_utils.get_resultant_wind(data["ua"], data["va"])

    # Step 2. Run Algorithm
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        if "time" in data.dims:
            data = data.squeeze("time")
        output = (
            jet_core_algorithms_components.run_jet_core_and_region_algorithm_on_one_day(
                data,
                jet_core_plev_limit,
                jet_core_ws_threshold,
                jet_boundary_ws_threshold,
                ws_drop_threshold,
                jet_core_lat_distance,
            )
        )
    else:
        output = data.groupby("time").map(
            jet_core_algorithms_components.run_jet_core_and_region_algorithm_on_one_day,
            (
                jet_core_plev_limit,
                jet_core_ws_threshold,
                jet_boundary_ws_threshold,
                ws_drop_threshold,
                jet_core_lat_distance,
            ),
        )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def penaortiz_et_al_2013(data):
    r"""
    This method follows a two step procedure for calculate local wind maxima days.
    This method returns 3 outputted variables:
        1. local_wind_maxima
        2. polar_front_jet
        3. subtropical_jet
    Each output is in a monthyear frequency.

    This method was first introduced in Pena-Ortiz et al. (2013) (https://doi.org/10.1002/jgrd.50305) and
    is described in section 2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing number of days per month with local wind maxima

    Notes
    -----
    Currently takes a long time i.e. 1.3 seconds per time unit with 8 plevs (i.e. 1.3 seconds per day)
    on AMD Ryzen 5 3600 6-core processor

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        uv_data = xr.open_dataset('path_to_uv_data')

        # Subset dataset to range used in original methodology (100-400 hPa)):
        uv_sub = uv_data.sel(plev=slice(100, 400))

        # Run algorithm:
        pena_outputs = jsmetrics.penaortiz_et_al_2013(uv_sub)

        # Produce a count of polar front jet values across all pressure levels
        pena_pfj_counts_all_levels = pena_outputs['polar_front_jet'].sum(('monthyear', 'plev'))

        # Use the polar front jet values as a mask to extract the jet windspeeds
        uv_sub["ws_by_monthyear"] = (
            data["ws"]
            .resample(time="MS")
            .sum()
            .rename({"time": "monthyear"})
        )
        pena_pfj_ws = pena_outputs.where(pena_outputs['polar_front_jet'] > 0)['ws_by_monthyear']

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
    r"""
    This method produces an event-based jet occurrence and jet center occurrence of JS.
    The outputs of this method will produce categorical values of three types:
        0. is not determined to be part of the jet
        1. is a jet occurence
        2. is jet core (upgraded from a jet occurence)

    This method was first introduced in Kuang et al (2014) (https://doi.org/10.1007/s00704-013-0994-x) and
    is described in section 2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.


    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    occurence_ws_threshold : int or float
        Threshold used to identify a jet-stream occurence point (default=30)

    Returns
    ----------
    output : xarray.Dataset
        Data containing jet-occurence and jet-centres (1 is occurence, 2 is core, 0 is no jet)

    Notes
    -----
    Currently takes a long time i.e. 2 seconds per time unit with 1 plev (i.e. 2 seconds per day) on AMD Ryzen 5 3600 6-core processor

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        uv_data = xr.open_dataset('path_to_uv_data')

        # Subset dataset to range used in original methodology (200 hPa)):
        uv_sub = uv_data.sel(plev=slice(200, 200))

        # Run algorithm:
        kuang_outputs = jsmetrics.kuang_et_al_2014(uv_sub, occurence_ws_threshold=30)

        # Extract only jet centers
        kuang_jet_centers = kuang_outputs.where(kuang_outputs['jet_ocurrence1_jet_centre2']==2)

        # Produce a count of jet centers across all pressure levels
        kuang_jet_centers_counts_all_levels = kuang_jet_centers.sum(('time', 'plev'))

        # Use the jet occurence values as a mask to extract the jet windspeeds
        kuang_jet_ws = kuang_outputs.where(kuang_outputs['jet_ocurrence1_jet_centre2'] > 0)['ws']

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
            if "time" in data.dims:
                data = data.squeeze("time")
            output = jet_core_algorithms_components.run_jet_occurence_and_centre_alg_on_one_day(
                data, occurence_ws_threshold
            )
        else:
            output = data.groupby("time").map(
                jet_core_algorithms_components.run_jet_occurence_and_centre_alg_on_one_day,
                (occurence_ws_threshold,),
            )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def jet_core_identification_algorithm(
    data, ws_core_threshold=40, ws_boundary_threshold=30
):
    r"""
    This method seperates jet cores based on boundary and windspeed threshold.
    Core are discovered where 8-cells are above boundary threshold

    This method is inspired by the method from Manney et al. (2011) (https://doi.org/10.5194/acp-11-6115-2011),
    which is described in Section 3.1 of that study.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    ws_core_threshold : int or float
        Threshold used for jet-stream core point (default=40)
    ws_boundary_threshold : int or float
        Threshold for jet-stream boundary point (default=30)

    Returns
    ----------
    output : xarray.Dataset
        Data containing the variable 'jet-core_id' (ID number relates to each unique core)

    Notes
    -----
    **Slow method**: currently takes a long time i.e. 7.6 seconds per time unit with 8 plevs
    (i.e. 7.6 seconds per day) on AMD Ryzen 5 3600 6-core processor

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        uv_data = xr.open_dataset('path_to_uv_data')

        # Subset dataset to range used in original methodology (100-400 hPa)):
        uv_sub = uv_data.sel(plev=slice(100, 400))

        # Run algorithm:
        jca_outputs = jsmetrics.jet_core_algorithms.jet_core_identification_algorithm(uv_sub, ws_core_threshold=40, ws_boundary_threshold=30)

    """
    if "plev" not in data.dims:
        data = data.expand_dims("plev")

    # Step 1. Run Jet-stream Core Idenfication Algorithm
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1:
        if "time" in data.dims:
            data = data.squeeze("time")
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
