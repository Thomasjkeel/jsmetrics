=====
Usage
=====

To use jsmetrics in a project::

.. code-block:: python

    import xarray as xr
    import jsmetrics

    # load windspeed data with u- and v- component wind.
    uv_data = xr.open_dataset(filename)

    # run Kuang et al. 2014 metric 
    k14 = jsmetrics.jetstream_algorithms.kuang_et_al_2014(uv_data)

.. image:: docs/_static/images/kuang_jet_centers.png
  :width: 360
  :align: center
  :alt: Kuang et al. 2014 Jet-core algorithm
