.. jsmetrics documentation master file, created by
   sphinx-quickstart on Wed Dec 29 17:38:41 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to jsmetrics's documentation!
==============================================

This is jsmetrics, a package containing implementations of various metrics and algorithms for identifying or characterising jet-streams
written in Python and built from xarray.
*WORK IN PROGRESS*

.. WRITE WHY JET-STREAM (maybe in blog, maybe in readme) -> heatwaves, beast from the east, climate proxy (put it all down)
.. At the foundation of studies that look at jet-streams is the metric used to describe or characterise it.

.. WRITE CURRENT PROGRESS WITH MODULDE in highlighted section near the top of this readme 

The philosophy of this package was to keep the methodology of each metric as close as possible to the given research paper's description of it (if not exact),
*but* to not limit the method to a given:
- time period,
- time unit (i.e. day, month, DJF),
- latitude/longitude resolution,
- latitude/longitude region (where possible),
- pressure level height,

All of these can be handled user-side.


.. 
        ALSO all algorithms have been broken down into various components and these components are not coupled to a given methodology.
        As such each can be used seperately and this allows users to rebuilt aspects of a methodology (e.g. to replace a filtering method)


Installation 
-------------
.. code-block:: bash
    
    pip install jsmetrics
    
Usage
-------------
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

DISCLAIMER
-------------
We have tried to replicate the various metrics based on the equations and details in the methodology as accurately as possible.
However, in some cases, we have chosen to exclude or alter parts of the methodology which reduce the resolution of the output (i.e. grouping into season or region) with the hope to preserve the parts of the method that specifically isolate a characteristics of the jet-stream at any inputted scale.
Again, any further subsetting is passed onto the user.
*If data input is at a daily resolution, part of the output should also be daily resolution.*  

Also note that, the data we used to test these metrics may have a different resolution to the one it was developed with.   

Finally, although these metric were found with a literature search, this is not an exaustive list of all methods used to identify or characterise the jet-stream or upper-level wind.
This project is very much a work in progress, so contributors are very welcome.

You can find details of each metric or algorithm here: `all metrics`_.

Where you can find my working-out (coming soon):
- I am hoping to make available all of my working out in jupyter-notebooks available soon (warning: these notebooks are not formatted) 
- I am also currently creating a verification notebook. 


Metrics & Algorithms
--------------------
See `all metrics`_ for specifications of each 'Complete' or 'In progress' metric and algorithm. For progress on their completion see `Status`_.


.. table::
   :align: left
   :widths: auto
   
   =============================================================================== ==============  ==  =============================================================================== ==============
   Metric/Algorithm                                                                `Status`_           Metric/Algorithm                                                                `Status`_                                                                                
   =============================================================================== ==============  ==  =============================================================================== ==============
   `Gallego et al. 2005 <http://link.springer.com/10.1007/s00382-005-0006-7>`_     To start            `Strong & Davis 2005 <http://doi.wiley.com/10.1029/2004GL022039>`_              To start
   `Koch et al. 2006 <https://onlinelibrary.wiley.com/doi/10.1002/joc.1255>`_      To verify           `Archer & Caldiera 2008 <http://doi.wiley.com/10.1029/2008GL033614>`_           To verify
   `Schiemann et al. 2009 <https://doi.org/10.1175/2008JCLI2625.1>`_               To verify           `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_   Complete
   `Manney et al. 2011 <https://acp.copernicus.org/articles/11/6115/2011/>`_       In progess*         `Allen et al. 2012 <http://www.nature.com/articles/nature11097>`_               To start
   `Barnes & Polvani 2013 <https://doi.org/10.1175/JCLI-D-12-00536.1>`_            To verify           `Pena-Ortiz et al. 2013 <http://doi.wiley.com/10.1002/jgrd.50305>`_             To verify      
   `Screen & Simmonds 2013 <http://doi.wiley.com/10.1002/grl.50174>`_              In progress*        `Kuang et al. 2014 <http://link.springer.com/10.1007/s00704-013-0994-x>`_       To verify            
   `Barnes & Polvani 2015 <https://doi.org/10.1175/JCLI-D-14-00589.1>`_            To verify           `Francis & Vavrus 2015 <https://doi.org/10.1088/1748-9326/10/1/014005>`_        Complete            
   `Cattiaux et al. 2016 <https://doi.wiley.com/10.1002/2016GL070309>`_            To verify           `Barnes & Simpson 2017 <https://doi.org/10.1175/JCLI-D-17-0299.1>`_             Complete            
   `Chenoli et al. 2017 <http://link.springer.com/10.1007/s00382-016-3102-y>`_     In progress         `Grise & Polvani 2017 <https://doi.org/10.1175/JCLI-D-16-0849.1>`_              Complete                        
   `Molnos et al. 2017  <https://doi.org/10.5194/esd-8-75-2017>`_                  In progress*        `Adam et al. 2018 <https://doi.org/10.5194/gmd-11-4339-2018>`_                  To start            
   `Bracegirdle et al. 2018 <https://doi.org/10.1175/JCLI-D-17-0320.1>`_           Complete            `Ceppi et al. 2018 <https://doi.org/10.1175/JCLI-D-17-0323.1>`_                 To verify            
   `Kern et al. 2018 <http://ieeexplore.ieee.org/document/8017585/>`_              To start*           `Rikus 2018 <http://dx.doi.org/10.1007/s00382-015-2560-y>`_                     In progress            
   `Kern & Westermann 2019 <https://doi.org/10.2312/vmv.20191321>`_                To start            `Kerr et al. 2020 <https://doi.org/10.1029/2020JD032735>`_                      To verify            
   `Maher et al. 2020 <https://doi.org/10.1007/s00382-019-05084-6>`_               To start            `Winters et al. 2020 <https://doi.org/10.1175/MWR-D-19-0353.1>`_                To start            
   `Martin 2021 <https://onlinelibrary.wiley.com/doi/10.1029/2020JD033668>`_       To start*           `Bosiger et al. 2022 <https://doi.org/10.5194/gmd-15-1079-2022>`_               To start            
   `Local Wave Activity <https://doi.org/10.1175/JAS-D-15-0194.1>`_                In progress*                        
   =============================================================================== ==============  ==  =============================================================================== ==============

* == help needed

.. _all metrics: https://github.com/Thomasjkeel/jsmetrics/blob/main/details_for_all_metrics.py
.. _Status: https://github.com/Thomasjkeel/jsmetrics/projects/1

.. 
        _also mention related references (i.e. Manney et al. )
        also Local Wave Activity (maybe martineu?)
        Gallego


.. Contributing
.. ------------
.. jsmetrics is in active development.

.. * If you're interested in participating in the development of jsmetrics by suggesting new features, new metrics or algorithms or report bugs, please leave us a message on the `issue tracker`_. There is also a chat room on gitter (|gitter|).

.. * If you would like to contribute code or documentation (which is greatly appreciated!), check out the `Contributing Guidelines`_ before you begin!

.. .. _issue tracker: https://github.com/Thomasjkeel/jsmetrics/issues
.. .. _Contributing Guidelines: https://github.com/Thomasjkeel/jsmetrics/blob/master/.github/CONTRIBUTING.rst


.. How to cite this package
.. ------------------------
.. If you wish to cite `jsmetrics` in a research publication, we kindly ask that you use the bibliographical reference information available through `Zenodo`


Credits
-------------

The layout and content of this project and was inspired by xclim (https://github.com/Ouranosinc/xclim) 
which contains other climate indices and metrics.

This package was created with Cookiecutter and the audreyr/cookiecutter-pypackage project template.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   metrics
   contributing
   authors
   history
   statement