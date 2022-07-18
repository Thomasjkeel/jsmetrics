==================
jsmetrics: Jet-stream metrics and algorithms
==================

|license| |pre-commit| |black| 

.. pypi| |conda| |coveralls| |codefactor|  |zenodo| |docs| 
----

This is jsmetrics, a library containing implementations of various metrics and algorithms for identifying or characterising jet-streams
written in Python and built from xarray.


.. [Table of Contents]
.. Disclaimer, table of metrics, Installation, Documentations, Contributing, How to cite, Project To-do's

DISCLAIMER
-------------
I have tried to replicate the various metrics based on the equations and details in the methodology as accurately as possible.

In some cases I have used a different dataset to the one used. 

This project is very much a work in progress, so contributors are very welcome. You  

Where you can find my working-out:
- I have included all of my working out in jupyter-notebooks available at: ... (warning: these notebooks have not been formatted nicely) 
- I am currently creating a verification notebook available at: ... where 


Table of the metrics
-------------
See [jsmetrics/details_for_all_metrics.py] for specifications of each 
For their progress see [Project 1]

.. table::
   :align: left
   :widths: auto
   
   ============================================================================== ==============  ==  ============================================================================== ==============
   Metric/Algorithm                                                               Status              Metric/Algorithm                                                               Status                                                                                
   ============================================================================== ==============  ==  ============================================================================== ==============
   `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_  Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_    Complete
   `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_  Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_    Complete
   `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_  Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_    Complete
   `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_  Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_    Complete
   `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_  Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_    Complete
   `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_  Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_    Complete
   `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_  Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_    Complete
   ============================================================================== ==============  ==  ============================================================================== ==============



Installation 
-------------
.. code-block:: python

    import jsmetrics
    import xarray as xr

    ua_data = xr.open_dataset(filename)
    w10 = jsmetrics.jetstream_metrics.woollings_et_al_2010(ua_data)

.. Documentation
.. -------------
.. The official documentation is at https://jsmetrics.readthedocs.io/

.. Contributing
.. ------------
.. jsmetrics is in active development and it's being used in production by climate services specialists.

.. * If you're interested in participating in the development of jsmetrics by suggesting new features, new indices or report bugs, please leave us a message on the `issue tracker`_. There is also a chat room on gitter (|gitter|).

.. * If you would like to contribute code or documentation (which is greatly appreciated!), check out the `Contributing Guidelines`_ before you begin!

.. .. _issue tracker: https://github.com/Thomasjkeel/jsmetrics/issues
.. .. _Contributing Guidelines: https://github.com/Thomasjkeel/jsmetrics/blob/master/.github/CONTRIBUTING.rst


.. How to cite this library
.. ------------------------
.. If you wish to cite `jsmetrics` in a research publication, we kindly ask that you use the bibliographical reference information available through `Zenodo`


Project To-Do's
-------------
- ADD: cf_xarray (see: https://cf-xarray.readthedocs.io/en/latest/index.html)
- ADD: pint (see: https://pint.readthedocs.io/en/stable/)
- ADD: DOI
- LOOK INTO timing/benchmarking the metrics (maybe in seperate github repo)
- TO SOLVE: dealing with data from different sources (some sort of data translator module or maybe included in tests)
  - for example what if 'v' or 'v-wind' is passed to func instead of 'va' (answer: cf-xarray)  
  - for example what if 'mbar' or 'model levels' instead of 'plev' (answer: pint)
- TO SOLVE: subsetting longitude if it wraps around 0-360
- CHECK: that methods using just U account for negative values (maybe need for abs() or not)

Credits
-------------

The layout and content of this project and was inspired by xclim (https://github.com/Ouranosinc/xclim) 
which contains other climate indices and metrics.

This package was created with Cookiecutter and the audreyr/cookiecutter-pypackage project template.

.. |license| image:: https://img.shields.io/badge/License-MIT-lightgray.svg?style=flt-square
        :target: https://github.com/Thomasjkeel/jsmetrics/blob/master/LICENSE
        :alt: License

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/python/black
        :alt: Python Black

.. |pre-commit| image:: https://results.pre-commit.ci/badge/github/Thomasjkeel/jsmetrics/master.svg
   :target: https://results.pre-commit.ci/latest/github/Thomasjkeel/jsmetrics/master
   :alt: pre-commit.ci status

.. .. |zenodo| image:: https://zenodo.org/badge/142608764.svg
..         :target: https://zenodo.org/badge/latestdoi/142608764
..         :alt: DOI

.. .. |docs| image:: https://readthedocs.org/projects/jsmetrics/badge
..         :target: https://jsmetrics.readthedocs.io/en/latest
..         :alt: Documentation Status

.. .. |pypi| image:: https://img.shields.io/pypi/v/jsmetrics.svg
..         :target: https://pypi.python.org/pypi/jsmetrics
..         :alt: Python Package Index Build

.. .. |conda| image:: https://img.shields.io/conda/vn/conda-forge/jsmetrics.svg
..         :target: https://anaconda.org/conda-forge/jsmetrics
..         :alt: Conda-forge Build Version
