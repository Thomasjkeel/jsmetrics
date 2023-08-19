===============
Examples of Use
===============

How to use jsmetrics in a project:


Using the jet statistics 
########################
...in combination with each other
---------------------------------


.. code-block:: python

    import jsmetrics
    import xarray as xr

    # Load in dataset with the variable 'ua' and coordinates: 'time', 'plev', 'lon' and 'lat':
    uv_data = xr.open_dataset('path_to_u_data')



Using the jet core algorithm outputs 
####################################

...as a spatial mask on other variables (such as windspeed)
-----------------------------------------------------------
Because all the jet core algorithm included in this package return 0 for regions not detected as the jet,
we can use xarray's `.where()` method to select a subset of another variable (i.e. windspeed)
within the boundaries of the detected jet.

.. code-block:: python

    import jsmetrics
    import xarray as xr

    # Load in dataset with the variables 'ua', 'va'; and coordinates: 'time', 'plev', 'lon' and 'lat'
    uv_data = xr.open_dataset('path_to_uv_data')

    # Subset dataset to range used in original methodology (100-500 hPa & 16.7-58.25 N, 42.5-220.5 E)):
    uv_sub = uv_data.sel(plev=slice(100, 500), lat=slice(16.7, 58.25), lon=slice(42.5, 220.5))

    # Run algorithm:
    schiemann_outputs = jsmetrics.jet_core_algorithms.schiemann_et_al_2009(uv_sub, ws_threshold=30)

    # Use the jet occurence values as a mask to extract the jet occurence windspeeds
    schiemann_jet_ws = schiemann.where(schiemann['jet_occurence'] > 0)['ws']


...to produce a count of jet cores:
------------------------------------
If you want to look at the frequency of jet locations and produce a map.

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
    schiemann_jet_counts_all_levels = schiemann['jet_occurence'].sum(('time', 'plev'))



