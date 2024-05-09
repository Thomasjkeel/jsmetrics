# -*- coding: utf-8 -*-

"""
    Statistics for isolating individual quantities synonymous with the jet stream from upper-level wind speed
    within a given time window (e.g. latitude, speed, width).

    The following statistics each return a xarray.Dataset and are ordered by paper publish year.
"""

# imports
import numpy as np
import xarray
from jsmetrics.utils import data_utils, spatial_utils, windspeed_utils
from jsmetrics.core.check_data import sort_xarray_data_coords
from . import jet_statistics_components

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


@sort_xarray_data_coords(coords=["lat", "lon"])
def archer_caldeira_2008(data):
    r"""
    This method calculates three monthly-averaged jet stream properties (windspeed, pressure and latitude) via integrated quantities from u-component wind.

    This method returns three properties:

    1. **weighted-average wind speed** -- jet stream wind speed (:math:`WS`), calculated by:

    .. math::
        WS_{i, j} =  \frac{\sum\limits_{k=400hPa}^{k=100hPa} m_{k} \times \sqrt{u^{2}_{i, j, k} + v^{2}_{i, j, k}}}{\sum\limits_{k=400hPa}^{k=100hPa} m_{k}}

    where :math:`u_{i,j,k}` and :math:`v_{i,j,k}` are the monthly-average horizontal wind components at grid point (i,j,k), and :math:`m_{k}` is the mass at level `k`.

    2. **mass flux weighted pressure** -- the average pressure of flows by the tropopause (:math:`P`), calculated by:

    .. math::
        P_{i, j} =  \frac{\sum\limits_{k=400hPa}^{k=100hPa} \left(m_{k} \times \sqrt{u^{2}_{i, j, k} + v^{2}_{i, j, k}}\right) \times p_k}{\sum\limits_{k=400hPa}^{k=100hPa} m_{k} \times \sqrt{u^{2}_{i, j, k} + v^{2}_{i, j, k}}}

    where :math:`p_k` is the pressure at level :math:`k`.

    3. **mass flux weighted latitude** -- Latitude of the Northern Hemisphere jet (:math:`L^{NH}`), calculated by:

    .. math::
        L_{i}^{NH} =  \frac{\sum\limits_{j=15°N}^{j=70°N} \left[\sum\limits_{k=400hPa}^{k=100hPa} \left(m_{k} \times \sqrt{u^{2}_{i, j, k} + v^{2}_{i, j, k}}\right) \right] \times \phi_{i,j}}{\sum\limits_{j=15N}^{j=70N} \sum\limits_{k=400hPa}^{k=100hPa} m_{k} \times \sqrt{u^{2}_{i, j, k} + v^{2}_{i, j, k}}}

    where :math:`\phi_{i,j}` is the grid cell latitude.

    This method was originally introduce in Archer & Caldiera (2008) (https://doi.org/10.1029/2008GL033614)
    and is described in Section 3 of that study.

    **Note:** this method does not explicitly limit inputted wind to 100-400 hPa, see 'Notes' for more information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing the three jet properties: 'mass_weighted_average_ws', 'mass_flux_weighted_pressure' and 'mass_flux_weighted_latitude'

    Notes
    -----
    While the initial methodology provides limits for pressure level (100-400 hPa), here the mass weighted outputs
    will be calculated for any pressure levels passed into the method.

    The latitude calculation is limited to 15-70N (as we only provide a way to extract the Northern Hemisphere jet),
    but you may find it easy enough to edit this method to calculate outputs for a different region.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u and v components:
        uv_data = xr.open_dataset('path_to_uv_data')

        # Subset dataset to range used in original methodology for the NH jet (100-400 hPa & 15-70 N)):
        uv_sub = uv_data.sel(plev=slice(100, 400), lat=slice(15, 70))

        # Run statistic:
        archer_outputs = jsmetrics.jet_statistics.archer_caldiera_2008(uv_sub)

        # Subset mass weighted wind by a windspeed threshold
        windspeed_threshold = 15 # remember this is for monthly averages
        strong_jet = archer_outputs['mass_weighted_average_ws'].where(lambda row: row > windspeed_threshold)

    """
    #  Step 1. Calculate monthly means
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

    #  Step 3. Calculate mass weighted average windspeed
    mass_weighted_average = jet_statistics_components.calc_mass_weighted_average(
        mon_mean, ws_col="ws"
    )

    #  Step 4. Calculate mass flux weighted pressure
    mass_flux_weighted_pressure = (
        jet_statistics_components.calc_mass_flux_weighted_pressure(
            mon_mean, ws_col="ws"
        )
    )

    #  Step 5. Calculate mass flux weighted latitude
    mass_flux_weighted_latitude = (
        jet_statistics_components.calc_mass_flux_weighted_latitude(
            mon_mean, lat_min=15, lat_max=75, ws_col="ws"
        )
    )

    #  Step 6. Assign the three new output variables to the original data
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
    r"""
    This method follows an in-text description of 4-steps describing the algorithm of jet-stream identification
    from Woollings et al. (2010).

    This method returns four outputs:
        1. **jet_lat** -- latitude of maximum speed within low-pass filtered zonally averaged wind profile
        2. **jet_speed** -- speed at the 'jet_lat'
        3. **ff_jet_lat** -- Fourier-filtered 'jet_lat' by season
        4. **ff_jet_speed** -- Fourier-filtered 'jet_speed' by season

    This method was first introduce in Woollings et al (2010) (http://dx.doi.org/10.1002/qj.625) and
    is described in section 2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package including Step 5 of the methodology.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    filter_freq : int
        number of days in filter (default=10 days)
    window_size : int
        number of days in window for Lancoz filter (default=61 days)

    Returns
    ----------
    output : xarray.Dataset
        Data containing the four output variables: 'jet_lat', 'jet_speed', 'ff_jet_lat' and 'ff_jet_speed'

    Notes
    -----
    In the original paper, a further step (Step 5) is carried out to express the values of jet latitude
    and jet speed anomalies from the seasonal cycle, this is shown in the 'Examples'

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to range used in original methodology (700-925 hPa & 20-70 N, 300-360 W)):
        ua_sub = ua.sel(plev=slice(700, 925), lat=slice(20, 70), lon=slice(300, 360))

        # Run statistic with a filter frequency and window size used in the original methodology:
        w10 = jsmetrics.jet_statistics.woollings_et_al_2010(ua_sub, filter_freq=10, window_size=61)

        # Express jet latitude and speed as anomalies from smoothed seasonal cycle (Step 5 of methodology)
        w10_seasonal_anomalies = w10.groupby('time.season').apply(lambda row: row['jet_lat'] - row['ff_jet_lat'])

    """
    # Initial translation of data from dataArray to dataset
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

    #  Step 6. Assign the new output variables to the original data
    time_dim = climatology["jet_speed"].dims[0]
    output = jet_statistics_components.assign_filtered_lats_and_ws_to_data(
        zonal_mean_lat_ws,
        fourier_filtered_lats.real,
        fourier_filtered_ws.real,
        dim=time_dim,
    )

    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def barnes_polvani_2013(data, filter_freq=10, window_size=41):
    r"""
    This method constructs the 'eddy-driven jet' by performing a pressure-weighted average of zonal winds. The winds
    are then low-pass frequency filtered at each grid point using a 10-day Lanczos filter with 41 weights by default.
    Finally a 0.01 degree quadratic function is fitted to the peak of the subsequent wind speed profile for each time step.

    This method returns three outputs:
        1. **jet_lat** -- latitude of maximum speed at an interval of 0.01 degree within the subseqeunt profile
        2. **jet_speed** -- speed at the 'jet_lat'
        3. **jet_width** -- full width at half of the maximum 'jet_speed' within the 0.01 quadratic fitted to the peak of the wind speed profile

    *Note:* 'jet_width' is undefined where the 'jet_speed' never drops below half maximum.

    This method was originally introduce in Barnes & Polvani (2013) https://doi.org/10.1175/JCLI-D-12-00536.1
    and is described in Section 2b and 2c of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    filter_freq : int
        number of days in filter (default=10 timeunits)
    window_size : int
        number of days in window for low-pass Lancoz filter (default=41 timeunits)

    Returns
    ----------
    output : xarray.Dataset
        Data containing the three output variables: 'jet_lat', 'jet_speed', 'jet_width'

    Notes
    -----
    Whereas the original analysis using this method focuses on three distinct sections of the globe,
    the method here does not make any distinction. Instead we highlight how to do subset the input data and
    calculate this metric for those three sections in 'Examples.

    This method was based on the method from Woollings et al. (2010) (http://dx.doi.org/10.1002/qj.625)

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to three sectors of the globe used in original methodology:
        ua_sh = ua.sel(plev=slice(700, 850), lat=slice(-90, 0)) # the Southern Hemisphere
        ua_na = ua.sel(plev=slice(700, 850), lat=slice(0, 90), lon=slice(300, 360) # the North Atlantic
        ua_np = ua.sel(plev=slice(700, 850), lat=slice(0, 90), lon=slice(135, 235) # the North Pacific

        # Run statistic with a filter frequency and window size used in the original methodology:
        bp13_sh = jsmetrics.jet_statistics.barnes_polvani_2013(ua_sh, filter_freq=10, window_size=41)
        bp13_na = jsmetrics.jet_statistics.barnes_polvani_2013(ua_na, filter_freq=10, window_size=41)
        bp13_np = jsmetrics.jet_statistics.barnes_polvani_2013(ua_np, filter_freq=10, window_size=41)
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

    #  Step 7. Assign the new output variables to the original data
    output = data.assign(
        {
            "jet_lat": (("time"), scaled_max_lats),
            "jet_speed": (("time"), scaled_max_ws),
            "jet_width": (("time"), jet_widths),
        }
    )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def grise_polvani_2014(data):
    r"""
    This method calculates the latitude of the midlatitude eddy-driven jet ('jet_lat') by finding the peak value of the input u-wind field.
    A polynomial fit is then applied to get an appropriate value of 'jet_lat' at a resolution 0.01 degrees.
    As opposed to the original method, this implementation also returns the speed at the 'jet_lat': the 'jet_speed'

    This method was originally introduce in Grise & Polvani (2014) https://doi.org/10.1002/2013GL058466
    and is described in Section 2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing the two outputs: 'jet_lat' and 'jet_speed'

    Notes
    -----
    This method was originally developed for the jet streams in the Southern Hemisphere.

    The original paper also includes two other metrics for zonal mean atmospheric circulation.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to range used in original methodology (850 hPa &  -65--30N)):
        ua_sub = ua.sel(plev=850, lat=slice(-65, -30))

        # Run statistic:
        gp16 = jsmetrics.jet_statistics.grise_polvani_2014(ua_sub)
    """
    if isinstance(data, xarray.DataArray):
        data = data.to_dataset()

    # Step 1. Expand time dimensions so we can map a function to the dataset properly
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1 and "time" not in data.dims:
        data = data.expand_dims("time")

    # Step 2. Calculate zonal-mean
    zonal_mean = windspeed_utils.get_zonal_mean(data)

    # Step 3. Get the 3 latitudes and speeds around max zonal wind-speed (e.g. lat-1, lat, lat+1)
    all_max_lats_and_ws = np.array(
        list(
            map(
                jet_statistics_components.get_3_latitudes_and_speed_around_max_ws,
                zonal_mean["ua"],
            )
        )
    )

    #  Step 4. Apply quadratic function to get max latitude at 0.01 degree resolution
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

    #  Step 5. Assign scaled max lats back to data
    output = data.assign(
        {
            "jet_lat": (("time"), scaled_max_lats),
            "jet_speed": (("time"), scaled_max_ws),
        }
    )
    return output


@sort_xarray_data_coords(coords=["lat", "lon"])
def barnes_polvani_2015(data):
    r"""
    This method calculates the jet positon and wind speed at that position by fitting a parabola around the
    maximum of zonally average u-component wind speed, using the magnitude at the maximum ('jet_speed') and latitude at that maximum ('jet_lat').

    This method was originally introduce in Barnes & Polvani (2015) https://doi.org/10.1175/JCLI-D-14-00589.1
    and is described in Section 2b of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing the two outputs: 'jet_lat' and 'jet_speed'

    Notes
    -----
    This methodology make an assumption that the a parabola can be fit to windspeed profile, so it performs quite different from
    other jet latitude methods available in the package in cases where the windspeed profile is more complex (multiple jets),
    and on data with finer temporal resolution in general.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to range used in original methodology (700-925 hPa &  30-70N,  130-10W)):
        ua_sub = ua.sel(plev=slice(700, 925), lat=slice(30, 70), lon=slice(230, 350))

        # Run statistic:
        bp15 = jsmetrics.jet_statistics.barnes_polvani_2015(ua_sub)
    """
    # Step 1. Get zonal mean
    zonal_mean = windspeed_utils.get_zonal_mean(data)

    # Step 2. Get jet lat and jet speed values and assign to output data
    if zonal_mean["time"].size == 1:
        if "time" in zonal_mean.dims:
            zonal_mean = zonal_mean.squeeze("time")
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
    r"""
    This method calculates two outputs: 'jet_lat' and 'jet_speed' which are defined as the latitude and speed of the 10-day-averaged
    maximum zonally-averaged wind speed.

    This method was originally introduce in Barnes & Simpson 2017 https://doi.org/10.1175/JCLI-D-17-0299.1
    and is described in Section 2b of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing the x outputs: 'jet_lat' and 'jet_speed'

    Notes
    -----
    The original methodology was intended to work on one pressure level (700 hPa) and on daily data, for the
    implementation included in this package, we have included methods to automatically average any inputted pressure levels
    and to return the data without 10-day averaging if data above the 10-day resolution is inputted.
    Instead, warnings are returned to user to help them use the method the way it was originally intended.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to range used in original methodology (700 hPa & North Atlantic & North Pacific)):
        ua_na = ua.sel(plev=700, lat=slice(0, 90), lon=slice(280, 350)) # North Atlantic
        ua_np = ua.sel(plev=700, lat=slice(0, 90), lon=slice(120, 230)) # North Pacific

        # Run statistic:
        bp17_na = jsmetrics.jet_statistics.barnes_simpson_2017(ua_na)
        bp17_np = jsmetrics.jet_statistics.barnes_simpson_2017(ua_np)

    """
    # Step 1. Run checks on the 'plev' coordinate. There should only be 1 for this method, but this implementation allows for multiple.
    if "plev" in data.dims:
        if data["plev"].count() == 1:
            data = data.isel(plev=0)
        else:
            print(
                "this metric was meant to only work on one plev, please subset plev to one value. For now taking the mean..."
            )
            data = data.mean("plev")
    data = data.mean("lon")

    # Step 1. Run check on time coordinate. The data should be able to be resampled into 10 days for this method.
    # Check 1. Is time coordinate is in data and is monotonically increasing.
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1 and "time" not in data.dims:
        data = data.expand_dims("time")
    if not data.indexes["time"].is_monotonic_increasing:
        raise IndexError("Data needs to have a montonic increasing index")
    # Check 2. Can that data be resampled into 10 days
    if not data["time"].size == 1:
        time_step_in_data = int((data["time"][1] - data["time"][0]).dt.days)
        if time_step_in_data <= 10:
            data = data.resample(time="10D").mean()
            time_step_in_data = 10
        else:
            print(
                f"Warning this method was developed for 10 day average and data has larger time-step than 10 days. Time step is {time_step_in_data} days"
            )
    data = data.dropna("time")  # Drop all NaN slices

    # Step 3. Calculate jet lat and jet speed.

    data = jet_statistics_components.calc_latitude_and_speed_where_max_ws(data)
    return data


@sort_xarray_data_coords(coords=["lat", "lon"])
def bracegirdle_et_al_2018(data):
    r"""
    This method calculates the seasonal and annual jet-stream position ('JPOS') and strength ('JSTR')
    by applying a 0.075 degrees cubic spline interpolation to a zonally-averaged wind climatology
    and selecting the maximum.

    This method was originally introduce in Bracegirdle et al (2018) https://doi.org/10.1175/JCLI-D-17-0320.1
    and is described in Section 2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing the four outputs: 'annual_JPOS', 'seasonal_JPOS', 'annual_JSTR' and 'seasonal_JSTR'

    Notes
    -----
    This method was originally developed for the jet streams in the Southern Hemisphere.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to range used in original methodology (850 hPa &  -75--10N)):
        ua_sub = ua.sel(plev=850, lat=slice(-75, -10))

        # Run statistic:
        b18 = jsmetrics.jet_statistics.bracegirdle_et_al_2018(ua_sub)
    """
    # Checks on plev
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

    # Step 1: Expand time dimensions so we can map a function to the dataset properly
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1 and "time" not in data.dims:
        data = data.expand_dims("time")

    #  Step 2. Calculate seasonal & annual climatologies
    seasonal_climatology = jet_statistics_components.get_climatology(data, "season")
    annual_climatology = jet_statistics_components.get_climatology(data, "year")

    #  Step 3. Get zonal mean from climatologies
    seasonal_zonal_mean = seasonal_climatology.mean("lon")
    annual_zonal_mean = annual_climatology.mean("lon")

    #  Step 4. Cubic spline interpolation to each climatology at latitude resolution of 0.075 degrees
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

    # Step 5. Assign jet-stream position (JPOS) and jet-stream strength (JSTR) back to data
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
    r"""
    This method calculates the jet latitude ('jet-lat') as defined by selecting the centroid of a zonally-averaged wind profile.

    The centroid is calculated by:

    .. math::
        \phi_{jet}  = \frac{\int_{30°}^{60°} \phi\bar{u}^2, d\phi}{\int_{30°}^{60°} \bar{u}^2, d\phi}

    This method has been slightly adapted to include a 'jet_speed' extraction (provided for this method
    in Screen et al. (2022) https://doi.org/10.1029/2022GL100523).

    **Note:** The implementation here does not explicit limit the centroid calculation to latitude between 20°-70°,
    instead this range is determined by the input data.

    This method was originally introduce in Ceppi et al (2018) https://doi.org/10.1175/JCLI-D-17-0323.1
    and is described in Section 2b of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    lon_resolution : numeric
        Resolution to use for longitude coord if size 1

    Returns
    ----------
    output : xarray.Dataset
        Data containing the three outputs: 'jet_lat', 'jet_speed', 'total_area_m2'

    Notes
    -----
    This method was improved by Zappa et al. (2018) https://doi.org/10.1029/2019GL083653,
    which includes exclusion of :math:`<0 m s^{-1}` u-wind. This method is also available in jsmetrics.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to range used in original methodology (850 hPa &  30-60N/S)):
        ua_na = ua.sel(plev=850, lat=slice(30, 60), lon=slice(300, 60)) # North Atlantic-European Sector
        ua_na = ua.sel(plev=850, lat=slice(30, 60), lon=slice(140, 240)) # North Pacific
        ua_sh = ua.sel(plev=850, lat=slice(-60, -30)) # Southern Hemisphere

        # Run statistic:
        c18_na = jsmetrics.jet_statistics.ceppi_et_al_2018(ua_na)
        c18_np = jsmetrics.jet_statistics.ceppi_et_al_2018(ua_np)
        c18_sh = jsmetrics.jet_statistics.ceppi_et_al_2018(ua_sh)
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

    # Step 4. Expand time dimension
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1 and "time" not in data.dims:
        data = data.expand_dims("time")
        zonal_mean = zonal_mean.expand_dims("time")

    # Step 5 (adapted from original methodology): Get nearest latitude actually in data to the one determined by metric
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

    # Step 6 (adapted from original methodology): Get speed of associated nearest latitude
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
    r"""
    This method calculates the jet latitude ('jet-lat') as defined by selecting the centroid of a zonally-averaged wind profile.

    The centroid is calculated by:

    .. math::
        \phi_{jet}  = \frac{\int_{20°}^{70°} \phi\bar{u}^2_0, d\phi}{\int_{20°}^{70°} \bar{u}^2_0, d\phi}

    .. math::
        u_0(\phi) = \max(0, u(\phi))

    **Note:** The implementation here does not explicit limit the centroid calculation to latitude between 20°-70°,
    instead this range is determined by the input data.

    This method has been slightly adapted to include a 'jet_speed' extraction (provided for this method
    in Screen et al. (2022) https://doi.org/10.1029/2022GL100523).

    This method was originally introduced in Zappa et al. 2018 https://doi.org/10.1029/2019GL083653
    and is described in Section 2.3 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    lon_resolution : numeric
        Resolution to use for longitude coord if size 1

    Returns
    ----------
    output : xarray.Dataset
        Data containing the three outputs: 'jet_lat', 'jet_speed', 'total_area_m2'

    Notes
    -----
    This method was adapted from and very similar to Ceppi et al (2018) https://doi.org/10.1175/JCLI-D-17-0323.1.

    This method is also used in Ayres & Screen, 2019 and Screen et al. 2022. Similar methods are used in:
    Chen et al. 2008; Ceppi et al. 2014, Ceppi et al. 2018.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to range used in original methodology (850 hPa &  20-70N, 140∘E-240∘E or 300-360 E))):
        ua_na = ua.sel(plev=850, lat=slice(20, 70), lon=slice(140, 240))
        ua_np = ua.sel(plev=850, lat=slice(20, 70), lon=slice(300, 360))


        # Run statistic:
        z18_na = jsmetrics.jet_statistics.zappa_et_al_2018(ua_na)
        z18_np = jsmetrics.jet_statistics.zappa_et_al_2018(ua_np)
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

    # Step 4. Expand time dimension
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    if data["time"].size == 1 and "time" not in data.dims:
        data = data.expand_dims("time")
        zonal_mean = zonal_mean.expand_dims("time")

    # Step 5 (adapted from original methodology): Get nearest latitude actually in data to the one determined by metric
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

    # Step 6 (adapted from original methodology): Get speed of associated nearest latitude
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
    r"""
    This method defines the latitude and speed of the jet-stream where the maximum zonal winds occur for
    each longitude and for each time unit (i.e. day). These values are then smoothed across the longitudes
    with a rectangular pulse (by default this has a width of 10 degrees).

    This method was originally introduced in Kerr et al. (2020) https://onlinelibrary.wiley.com/doi/10.1029/2020JD032735
    and is described in Section 2.4.2 of that study.

    Please see 'Notes' below for any additional information about the implementation of this method
    to this package.

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua', and the coordinates: 'lon', 'lat', 'plev' and 'time'.

    Returns
    ----------
    output : xarray.Dataset
        Data containing the two outputs: 'jet_lat' and 'smoothed_jet_lat'

    Notes
    -----
    This method was based on the method from Barnes and Fiore (2013) https://doi.org/10.1002/grl.50411

    The implementation here returns both smoothed and unsmoothed jet latitude outputs.

    Examples
    --------
    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Load in dataset with u component wind:
        ua_data = xr.open_dataset('path_to_u_data')

        # Subset dataset to range used in original methodology (700 hPa &  20-70N,  300-360W)):
        ua_sub = ua.sel(plev=700, lat=slice(20, 70), lon=slice(300, 360))

        # Run statistic:
        k20 = jsmetrics.jet_statistics.kerr_et_al_2020(ua_sub)
    """
    # Checks on plev coordinate
    if "plev" in data.dims:
        if data["plev"].count() == 1:
            data = data.isel(plev=0)
        else:
            print(
                "this metric was meant to only work on one plev, please subset plev to one value. For now taking the mean..."
            )
            data = data.mean("plev")

    # Step 1. Calculateed smoothed jet lats by lon
    if "time" not in data.coords:
        raise KeyError("Please provide a time coordinate for data to run this metric")
    elif data["time"].size == 1:
        if "time" in data.dims:
            data = data.squeeze("time")
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
