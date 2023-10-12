=======================
Statistics & Algorithms
=======================
:code:`jsmetrics` contains metrics of three types: 
   1. **jet statistics** -- methods that extract a single quantity describing the jet stream (i.e. latitude, speed, width)
   2. **jet core algorithms** -- methods that return a multi-dimensional mask at coordinates identified to be the jet streams.
   3. **waviness metrics** -- methods that characterise the sinousity of the upper-level winds (acts as a proxy for jet stream waviness).


For specification details of each metric's DOI and original study area please see `all metrics`_.

For progress on their completion see `issues`_.

.. _all metrics: https://github.com/Thomasjkeel/jsmetrics/blob/main/jsmetrics/details_for_all_metrics.py
.. _issues: https://github.com/Thomasjkeel/jsmetrics/issues

Jet statistics
##############
Metrics and indices used to characterise features of jet-streams (e.g. location, speed, waviness and width).
All methods are implementation from existing research which have been translated to Python for this package.

**Overview table:**

.. table::
   :align: left
   :widths: auto
   
   
   =============================================================================== ==============  ==  =============================================================================== ==============
   Metric                                                                          Status              Metric                                                                          Status                                                                                
   =============================================================================== ==============  ==  =============================================================================== ==============
   `Archer & Caldiera 2008 <http://doi.wiley.com/10.1029/2008GL033614>`_           Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_   Complete
   `Barnes & Polvani 2013 <https://doi.org/10.1175/JCLI-D-12-00536.1>`_            Complete            `Screen & Simmonds 2013 <http://doi.wiley.com/10.1002/grl.50174>`_              In progress*    
   `Grise & Polvani 2014 <https://doi.org/10.1002/2013GL058466>`_                  Complete            `Barnes & Polvani 2015 <https://doi.org/10.1175/JCLI-D-14-00589.1>`_            Complete
   `Barnes & Simpson 2017 <https://doi.org/10.1175/JCLI-D-17-0299.1>`_             Complete            `Chenoli et al. 2017 <http://link.springer.com/10.1007/s00382-016-3102-y>`_     In progress
   `Molnos et al. 2017  <https://doi.org/10.5194/esd-8-75-2017>`_                  In progress*        `Bracegirdle et al. 2018 <https://doi.org/10.1175/JCLI-D-17-0320.1>`_           Complete
   `Ceppi et al. 2018 <https://doi.org/10.1175/JCLI-D-17-0323.1>`_                 Complete            `Zappa et al. 2018 <https://doi.org/10.1029/2019GL083653>`_                     Complete
   `Rikus 2018 <http://dx.doi.org/10.1007/s00382-015-2560-y>`_                     In progress         `Kerr et al. 2020 <https://doi.org/10.1029/2020JD032735>`_                      To verify
   =============================================================================== ==============  ==  =============================================================================== ==============

\* == help needed

**Documentation:**

.. automodule:: jsmetrics.metrics.jet_statistics
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Jet core algorithms
###################
Methods that isolate parts of the atmospheric column associated with jet-streams based on various techniques from the literature.
In each case, these algorithms output a mask of the input data, and so can be used to further characterise variables at the same location
as the jet core mask. Algorithms are treated seperately to jet statistcs, as in general, the jet statistic are used to summarise information
about the jet-stream, the jet core algorithms simply identify it. All methods are implementation from existing research which have been
translated to Python for this package.


**Overview table:**

.. table::
   :align: left
   :widths: auto
   
   =============================================================================== ==============  ==  =============================================================================== ==============
   Algorithm                                                                       Status              Algorithm                                                                       Status                                                                                
   =============================================================================== ==============  ==  =============================================================================== ==============
   `Koch et al. 2006 <https://onlinelibrary.wiley.com/doi/10.1002/joc.1255>`_      Complete            `Schiemann et al. 2009 <https://doi.org/10.1175/2008JCLI2625.1>`_               Complete          
   `Manney et al. 2011 <https://acp.copernicus.org/articles/11/6115/2011/>`_       To verify           `Pena-Ortiz et al. 2013 <http://doi.wiley.com/10.1002/jgrd.50305>`_             To verify      
   `Kuang et al. 2014 <http://link.springer.com/10.1007/s00704-013-0994-x>`_       To verify                 
   =============================================================================== ==============  ==  =============================================================================== ==============

\* == help needed

**Documentation:**

.. automodule:: jsmetrics.metrics.jet_core_algorithms
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Waviness metrics
################
Methods that characterise 'waviness' or 'sinousity' in parts of the atmospheric column based on various techniques from the literature.
All methods are implementation from existing research which have been translated to Python for this package.


**Overview table:**

.. table::
   :align: left
   :widths: auto
   
   =============================================================================== ==============  ==  =============================================================================== ==============
   Algorithm                                                                       Status              Algorithm                                                                       Status                                                                                
   =============================================================================== ==============  ==  =============================================================================== ==============
   `Francis & Vavrus 2015 <https://doi.org/10.1088/1748-9326/10/1/014005>`_        Complete            `Cattiaux et al. 2009 <https://doi.org/10.1002/2016GL070309>`_                  To verify          
   `Local Wave Activity <https://doi.org/10.1175/JAS-D-15-0194.1>`_                In progress*
   =============================================================================== ==============  ==  =============================================================================== ==============

\* == help needed

**Documentation:**

.. automodule:: jsmetrics.metrics.waviness_metrics
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :noindex:

See the :ref:`Metric sub-components` for more detail about the implementation of the metrics included in this package.