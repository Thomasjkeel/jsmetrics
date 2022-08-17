.. jsmetrics documentation master file, created by
   sphinx-quickstart on Wed Dec 29 17:38:41 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to jsmetrics's documentation!
==============================================

This is jsmetrics, a library containing implementations of various metrics and algorithms for identifying or characterising jet-streams
written in Python and built from xarray.


WRITE WHY JET-STREAM (maybe in blog, maybe in readme) -> heatwaves, beast from the east, climate proxy (put it all down)

WRITE CURRENT PROGRESS WITH MODULDE in highlighted section near the top of this readme 

NOT DEPENDENT ON time, lat/lon resolution, lat/lon bbox, plev height (all changable)
ALSO can use components to rebuilt aspects of the methodology 
.. [Table of Contents]
.. Disclaimer, table of metrics, Installation, Documentations, Contributing, How to cite, Project To-do's

DISCLAIMER
-------------
I have tried to replicate the various metrics based on the equations and details in the methodology as accurately as possible.

In some cases I have used a different dataset to the one used. 

This project is very much a work in progress, so contributors are very welcome. You  

Details provided in: `all metrics`_ is not exact as in some cases ... Most algorithms can be used at different pressure-levels etc.  

Where you can find my working-out:
- I have included all of my working out in jupyter-notebooks available at: ... (warning: these notebooks have not been formatted nicely) 
- I am currently creating a verification notebook available at: ... where 


Table of the metrics
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
   `Cattiaux et al. 2016 <https://doi.wiley.com/10.1002/2016GL070309>`_            To verify           `Barnes & Simpson 2017 <https://doi.org/10.1175/JCLI-D-17-0299.1>`_             To verify            
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

.. _all metrics: https://github.com/Thomasjkeel/jsmetrics/blob/write-docs/jsmetrics/details_for_all_metrics.py
.. _Status: https://github.com/Thomasjkeel/jsmetrics/projects/1


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   installation
   usage
   metrics
   contributing
   authors
   history

