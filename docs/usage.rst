===============
Examples of Use
===============

:code:`jsmetrics` is designed to be easy to use and should integrate seemlessly with `*xarray* <https://docs.xarray.dev/en/stable/>`_ 
(which is a Python library useful for dealing with earth science datasets).
An extensive knowledge of Python or *xarray* is **not** required to use *jsmetrics*, although it will help you use the package
more effectively if you wish to run some of the more advanced use cases. 

.. note:: 
    To run any metric in :code:`jsmetrics` the syntax will be something like:

    .. code-block:: python

        import jsmetrics
        import xarray as xr

        # Use xarray to load in NetCDF or GRIB format data
        your_data = xr.open_dataset('path_to_your_data.nc')

        # Run a metric on your data and store the outputs 
        output = jsmetrics.<jet_module>.<jet_metric>(your_data)


*jsmetrics* provides three :ref:`types <Statistics & Algorithms>` of metric, we provide examples for each of them:
    1. :ref:`Jet statistics <1. Using the jet statistics>` 
    2. :ref:`Jet core algorithms <2. Using the jet core algorithms>`
    3. :ref:`Waviness metrics <3. Using the waviness metrics>`

*Please note that we also provide some examples in a jupyter notebook format available* `here <https://github.com/Thomasjkeel/jsmetrics-examples>`_.

1. Using the jet statistics 
###########################
...in combination with each other
---------------------------------


.. code-block:: python

    import jsmetrics
    import xarray as xr

    # Load in dataset with the variable 'ua' and coordinates: 'time', 'plev', 'lon' and 'lat':
    uv_data = xr.open_dataset('path_to_u_data')



2. Using the jet core algorithms 
################################

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


3. Using the waviness metrics 
#############################

