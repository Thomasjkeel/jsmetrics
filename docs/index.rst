.. jsmetrics documentation master file, created by
   sphinx-quickstart on Wed Dec 29 17:38:41 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to jsmetrics' documentation!
==============================================
*jsmetrics* is an open-source Python package containing implementations of various statistics and algorithms which have been
developed to identify or characterise :ref:`jet streams<What are jet streams?>`.

The package is built on top of `*xarray* <https://docs.xarray.dev/en/stable/>`_ and currently contains 17 methods,
consisting of jet statistics, waviness metrics and jet core algorithms (described :ref:`here <Statistics & Algorithms>`).
As this is an ongoing project, we are always in the process of finding and implementing new methods.
You can find a full list of methods and their current progress state :ref:`here <Jet stream metrics>`

.. note::
  - preprint now available on `EGUsphere <https://egusphere.copernicus.org/preprints/2023/egusphere-2023-661/>`_
  - example jupyter notebooks on `Github <https://github.com/Thomasjkeel/jsmetrics-examples>`_


Jet stream metrics
------------------
.. table::
   :align: left
   :widths: auto

   =============================================================================== ==============  ==  =============================================================================== ==============
   Statistic/Algorithm                                                             `Status`_           Statistic/Algorithm                                                             `Status`_
   =============================================================================== ==============  ==  =============================================================================== ==============
   `Gallego et al. 2005 <http://link.springer.com/10.1007/s00382-005-0006-7>`_     To start            `Strong & Davis 2005 <http://doi.wiley.com/10.1029/2004GL022039>`_              To start
   `Koch et al. 2006 <https://onlinelibrary.wiley.com/doi/10.1002/joc.1255>`_      Complete            `Archer & Caldiera 2008 <http://doi.wiley.com/10.1029/2008GL033614>`_           Complete
   `Schiemann et al. 2009 <https://doi.org/10.1175/2008JCLI2625.1>`_               Complete            `Woollings et al. 2010 <https://onlinelibrary.wiley.com/doi/10.1002/qj.625>`_   Complete
   `Manney et al. 2011 <https://acp.copernicus.org/articles/11/6115/2011/>`_       To verify           `Allen et al. 2012 <http://www.nature.com/articles/nature11097>`_               To start
   `Barnes & Polvani 2013 <https://doi.org/10.1175/JCLI-D-12-00536.1>`_            Complete            `Pena-Ortiz et al. 2013 <http://doi.wiley.com/10.1002/jgrd.50305>`_             To verify
   `Screen & Simmonds 2013 <http://doi.wiley.com/10.1002/grl.50174>`_              In progress*        `Grise & Polvani 2014 <https://doi.org/10.1002/2013GL058466>`_                  Complete
   `Kuang et al. 2014 <http://link.springer.com/10.1007/s00704-013-0994-x>`_       To verify           `Barnes & Polvani 2015 <https://doi.org/10.1175/JCLI-D-14-00589.1>`_            Complete
   `Francis & Vavrus 2015 <https://doi.org/10.1088/1748-9326/10/1/014005>`_        Complete            `Cattiaux et al. 2016 <https://doi.wiley.com/10.1002/2016GL070309>`_            To verify
   `Barnes & Simpson 2017 <https://doi.org/10.1175/JCLI-D-17-0299.1>`_             Complete            `Chenoli et al. 2017 <http://link.springer.com/10.1007/s00382-016-3102-y>`_     In progress
   `Bracegirdle et al. 2018 <https://doi.org/10.1175/JCLI-D-17-0320.1>`_           Complete            `Ceppi et al. 2018 <https://doi.org/10.1175/JCLI-D-17-0323.1>`_                 Complete
   `Rikus 2018 <http://dx.doi.org/10.1007/s00382-015-2560-y>`_                     In progress         `Zappa et al. 2018 <https://doi.org/10.1029/2019GL083653>`_                     Complete
   `Kerr et al. 2020 <https://doi.org/10.1029/2020JD032735>`_                      To verify           `Maher et al. 2020 <https://doi.org/10.1007/s00382-019-05084-6>`_               To start
   `Bosiger et al. 2022 <https://doi.org/10.5194/gmd-15-1079-2022>`_               To start            `Local Wave Activity <https://doi.org/10.1175/JAS-D-15-0194.1>`_                In progress*
   =============================================================================== ==============  ==  =============================================================================== ==============

\* == help needed

For details of each 'Complete' or 'In progress' metric and algorithm, see our `specification file`_.
We track the progress of each method using a GitHub Kanban, and you can see their `status`_.

.. _specification file: https://github.com/Thomasjkeel/jsmetrics/blob/main/jsmetrics/details_for_all_metrics.py
.. _status: https://github.com/Thomasjkeel/jsmetrics/projects/1


.. Examples
.. -------------
.. For examples please check out the :ref:`Examples of Use`.

.. Some example jupyter notebooks are also made available `here <https://github.com/Thomasjkeel/jsmetrics-examples>`_

.. Disclaimer
.. ----------
.. We have tried to replicate the various metrics based on the equations and details in the methodology as accurately as possible.
.. However, in some cases, we have chosen to exclude or alter parts of the methodology which reduce the resolution of the output (i.e. grouping into season or region) with the hope to preserve the parts of the method that specifically isolate a characteristics of the jet-stream at any inputted scale.
.. Again, any further subsetting is passed onto the user.
.. *If data input is at a daily resolution, part of the output should also be daily resolution.*

.. Also note that, the data we used to test these metrics may have a different resolution to the one it was developed with.

.. Finally, although these metric were found with a literature search, this is not an exaustive list of all methods used to identify or characterise the jet-stream or upper-level wind.
.. This project is very much a work in progress, so contributors are very welcome.

.. You can find details of each metric or algorithm here: `all metrics`_.

..
        _also mention related references (i.e. Manney et al. )
        also Local Wave Activity (maybe martineu?)
        Gallego

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


.. Credits
.. -------------

.. The layout and content of this project and was inspired by xclim (https://github.com/Ouranosinc/xclim)
.. which contains other climate indices and metrics.

.. This package was created with Cookiecutter and the audreyr/cookiecutter-pypackage project template.

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:

   installation
   statement
   usage
   metrics
   components
   contributing
   authors
   history
