# -*- coding: utf-8 -*-

"""
    Methods that return a mask of coordinates related to the jet location, e.g., identifying the maximum
    wind speed throughout the horizontal and/or vertical plane within a given time window.

    The following algorithms each return a xarray.Dataset and are ordered by paper publish year.
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
    This method follows a two-step procedure used to detect jet-event occurences (here: 'jet_events_ws').

    The weighted average windspeed (:math:`\alpha vel`) for the jet events is calculated as follows:

    .. math::
        \alpha vel =  \frac{1}{p2-p1} \int_{p1}^{p2} (u^2+v^2)^{1/2} \,dp

    where :math:`p1`, :math:`p2` is min, max pressure level.

    After calculating :math:`\alpha vel`, in a second step a windspeed threshold to isolate jet events (the default is :math:`30 m s^{-1}`).

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
    output : xr.Dataset
        Data containing the output variable 'jet_events_ws'

    Notes
    -----
    The original methodology uses windspeed between 100-400 hPa to calculated the weighted average
    and 30 meters per second as the windspeed threshold.

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
    # Step 1. Check there are at least two plev coordinates for weighted average.
    if data["plev"].count() < 2:
        raise ValueError(
            "Need at least two plevs to calculate weighted average windspeed"
        )

    # Step 2: get all pressure levels (hPa) as list
    all_plevs_hPa = jet_core_algorithms_components.get_all_hPa_list(data)

    # Step 3: get weighted sum windspeed
    sum_weighted_ws = jet_core_algorithms_components.get_sum_weighted_ws(
        data, all_plevs_hPa
    )

    # Step 4: calculate average weighted
    weighted_average_ws = jet_core_algorithms_components.get_weighted_average_ws(
        sum_weighted_ws, all_plevs_hPa
    )

    # Step 5: Apply windspeed threshold to get jet event array
    jet_events = weighted_average_ws.where(weighted_average_ws >= ws_threshold)
    jet_events = jet_events.fillna(0.0)

    # Step 6: Turn jet event array into dataset
    jet_event_ds = jet_events.rename("jet_events_ws").to_dataset()
    return jet_event_ds


@sort_xarray_data_coords(coords=["lat", "lon"])
def schiemann_et_al_2009(data, ws_threshold=30, u_threshold=0):
    r"""
    This method detects 'jet occurrences', whereby each jet occurence is detected based
    on three rules applied to inputted wind speed (:math:`V = [u, v]`):

    1. :math:`|V|` is a local maxima in latitude and altitude plane
    2. :math:`|V| \ge 30 m s^{-1}`
    3. :math:`u \ge 0 m s^{-1}`.

    The implementation of this method here allows you to edit the :math:`|V|` threshold
    (by changing 'ws_threshold'), and :math:`u` threshold (by changing 'u_threshold').

    This method was originally introduce in Schiemann et al 2009 (https://doi.org/10.1175/2008JCLI2625.1)
    and is described in Section 2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    ws_threshold : int or float
        Windspeed threshold used to extract jet maxima from resultant windspeed (default: 30 m/s)
    u_threshold : int or float
        Windspeed threshold used to extract u-component wind speed (default: 0 m/s)

    Returns
    ----------
    output : xr.Dataset
        Data containing the two output variables: 'ws' and 'jet_occurence'

    Notes
    -----
    While the original method is built on a four dimension slice of wind speed (time, lat, lon, plev),
    This implementation will work where there is only one pressure level, so a 3-d slice (time, lat, lon).

    **Slow method:** due to the nature of this method, it currently takes a moderately long time to run,
    i.e. 7.6 seconds per time unit on AMD Ryzen 5 3600 6-core processor.

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

        # Get jet occurence for first day in data
        schiemann_jet_occurence_first_day = schiemann_outputs['jet_occurence'].isel(time=0).max('plev')

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
        output = jet_core_algorithms_components.get_local_jet_occurence(
            data, ws_threshold=ws_threshold, u_threshold=u_threshold
        )
    else:
        output = data.groupby("time").map(
            jet_core_algorithms_components.get_local_jet_occurence,
            (
                ws_threshold,
                u_threshold,
            ),
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
    check_diagonals=False,
):
    r"""
    This method detects jet cores (within an altitude range see 'jet_core_plev_limit') and a boundary region surrounding
    those cores based on two windspeed thresholds. Two checks are applied after initial detection of cores to check whether
    boundaries with more then one core are part of the same feature (the default threshold for these boundaries is 30 m/s,
    see 'jet_boundary_ws_threshold'). The two checks seperate cores based on whether the cores are more than a certain distance apart
    (default is 15 degrees, see 'jet_core_lat_distance') and whether the windspeed between two given cores does not drop
    below a windspeed threshold (default is 25 m/s, see 'ws_drop_threshold')

    This method returns four outputs
        1. **jet_core_mask** -- Regions within each latitude-altitude slice that are local maxima and have windspeeds above the 'jet_core_ws_threshold'
        2. **jet_region_mask** -- Regions above, below, left and right of any given jet core with windspeed above the 'jet_boundary_ws_threshold'
        3. **jet_region_contour_mask** -- All contigious regions of windspeeds emcompassing a jet core above the 'jet_boundary_ws_threshold' (i.e. not just above, below, left and right)
        4. **ws** -- Resultant wind speed calculated from 'ua', 'va' inputs.

    This method was originally introduce in Manney et al. (2011) (https://doi.org/10.5194/acp-11-6115-2011),
    and is described in Section 3.1 of that study. This method is also known as the JETPAC (Jet and
    Tropopause Products for Analysis and Characterization) software package, and available in its original
    form at request to NASA JPL.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    jet_core_plev_limit: tuple or array
        Sequence of two values relating to the pressure level limit of the jet cores (original paper uses 100-400 hPa)
    jet_core_ws_threshold : int or float
        Threshold used for jet cores (default=40 m/s)
    jet_boundary_ws_threshold : int or float
        Threshold for jet boundaries (default=30 m/s)
    ws_drop_threshold : int or float
        Threshold for drop in windspeed along direct interpolated path between cores (default: 25 m/s)
    jet_core_lat_distance : int or float
        Threshold for maximum distance between cores to be counted the same (default: 15 degrees)
    check_diagonals : bool
        Whether to check the diagonal edges of each maxima in the latitude-altitude plane. The original method does not include this (default: False).

    Returns
    ----------
    output : xarray.Dataset
        Data containing the four output variables: 'ws', 'jet_region_mask', 'jet_region_contour_mask', and 'jet_core_mask'

    Notes
    -----
    The implementation of this method varies slightly from the original, because this method will return a mask rather
    than dynamical values. The intention of this, is to allow these masks to be used to subset other variables such as windspeed
    (see 'Examples' for demonstration of how to use these masks).

    There is an update to this method introduced in Manney & Hegglin 2018 to include physically-based method to extract the
    subtropical jet is identified (and thus distinguished from polar jets). This update is not included in jsmetrics.

    'jet_region_above_ws_threshold_mask' is provided here, but not explicitly in the original methodology. It is meant as
    an alternative to using a contour to check which regions encompass jet cores.

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
    # Step 1. Check plev and time coordinate in data
    if "plev" not in data.dims:
        data = data.expand_dims("plev")
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")

    # Step 2. Check a pressure level limit is provided by the user
    if not jet_core_plev_limit:
        raise KeyError(
            "Please provide a pressure level limit for jet cores returned by this metric. As an example the original methodology used a limit of 100-400 hPa. To replicate this, pass the parameter jet_core_plev_limit=(100, 400)."
        )

    # Step 3. Calculate wind speed from ua and va components.
    data["ws"] = windspeed_utils.get_resultant_wind(data["ua"], data["va"])

    # Step 4. Run Algorithm and return outputs
    if data["time"].size == 1:
        output = (
            jet_core_algorithms_components.run_jet_core_and_region_algorithm_on_one_day(
                data,
                jet_core_plev_limit,
                jet_core_ws_threshold,
                jet_boundary_ws_threshold,
                ws_drop_threshold,
                jet_core_lat_distance,
                check_diagonals,
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
                check_diagonals,
            ),
        )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def penaortiz_et_al_2013(data, ws_threshold=30):
    r"""
    This method follows a two step procedure for calculating local wind maxima and then subcategorising the local maxima into
    two distinct jet masks: the Subtropical Jet (STJ) and Polar Front Jet (PFJ).

    This method returns 4 outputs:
        1. **local_wind_maxima** -- Binary mask of local wind maxima
        2. **local_wind_maxima_by_monthyear** -- Same as above, but for monthyear frequency
        3. **polar_front_jet** -- Binary mask of the PFJ (by monthyear).
        4. **subtropical_jet** -- Binary ask of the STJ (by monthyear).

    This method was first introduced in Pena-Ortiz et al. (2013) (https://doi.org/10.1002/jgrd.50305) and
    is described in section 2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    ws_threshold : int or float
        windspeed threshold to apply (default=30 ms-1)

    Returns
    ----------
    output : xarray.Dataset
        Data containing the four output variables: 'local_wind_maxima', 'local_wind_maxima_by_monthyear', 'polar_front_jet', and 'subtropical_jet'

    Notes
    -----
    See Table 1 in the respective paper for the categories used to seperate the STJ and PFJ.
    The STJ is only seperated in DJF for the Northern Hemisphere.

    **Slow method**: currently takes a long time i.e. 1.3 seconds per time unit with 8 plevs (i.e. 1.3 seconds per day)
    on a AMD Ryzen 5 3600 6-core processor.

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
        jet_core_algorithms_components.get_local_wind_maxima_by_timeunit,
        (ws_threshold,),
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
    This method produces an event-based jet occurrences and jet center occurrences of the jet stream
    in a given atmospheric column.

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
    **Slow method**: currently takes a long time i.e. 2 seconds per time unit with 1 plev (i.e. 2 seconds per day) on AMD Ryzen 5 3600 6-core processor

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
    # Step 1. Check plev and time coordinate in data
    if "plev" not in data.dims:
        data = data.expand_dims("plev")
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")

    # Step 2. Calculate wind speed from ua and va components.
    data["ws"] = windspeed_utils.get_resultant_wind(data["ua"], data["va"])

    # Step 2. Run Jet-stream Occurence and Centre Algorithm and return outputs
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
    This method extract seperate jet cores based on boundary and core windspeed thresholds.

    The output variable of this method includes two types:
        - 0 -- regions not determined to be part of the jet
        - 1-n -- Seperate jet core regions seperated by one condition: if two cores in the same region are more than 15 degrees of latitude away

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

        # Subset dataset to range appropriate for methodology (100-400 hPa)):
        uv_sub = uv_data.sel(plev=slice(100, 400))

        # Run algorithm:
        jca_outputs = jsmetrics.jet_core_algorithms.jet_core_identification_algorithm(uv_sub, ws_core_threshold=40, ws_boundary_threshold=30)

    """
    # Step 1. Check plev coordinate in data
    if "plev" not in data.dims:
        data = data.expand_dims("plev")

    # Step 2. Run Jet-stream Core Idenfication Algorithm and return outputs
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
