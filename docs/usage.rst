===============
Examples of Use
===============

To use jsmetrics in a project::

.. code-block:: python

 import xarray as xr
 import jsmetrics

 # load windspeed data with u- and v- component wind.
 uv_data = xr.open_dataset(filename)

 # run Woollings et al. 2010 metric
 w10 = jsmetrics.metrics.jet_statistics.woollings_et_al_2010(uv_data)

 print(w10['jet_lat'])
 print(w10['jet_speed'])

 # run Kuang et al. 2014 metric. NOTE: may take a long time after you have more than 50 time steps.
 k14 = jsmetrics.metrics.jet_core_algorithms.kuang_et_al_2014(uv_data)
 print(k14['jet_center'].sel(time=0))


.. image:: docs/_static/images/all_metrics_jetlat_circbar_w_errorbars.png
  :width: 560
  :align: center
  :alt: Jet latitude circbars with errorbars

*Estimation of North Pacific mean jet latitude by month with 1-stdev errorbars. Data is monthly ERA5 700-850 hPa u-wind between 1980-2020.*

.. image:: docs/_static/images/jet_core_algorithm_comparions_NA_5_texas2021.png
  :width: 560
  :align: center
  :alt: Comparison of jet core algorithms during Feb 2021 Texas Cold Wave

*Comparison of jet core algorithms estimation of the 6-hourly jet position. Data is 6-hourly ERA5 100-500 hPa u-v-wind.*


.. image:: docs/_static/images/all_jet_lats_stj_pfj_npac_maps_more_metrics.png
  :width: 560
  :align: center
  :alt: STJ and PFJ by metric and longitude

*By latitude estimation of the jet latitude of the subtropical and polar jet stream. Data is monthly ERA5 differenced-250 hPa (orange) and 700-850 hPa (blue) u-wind between 1980-2020.*

