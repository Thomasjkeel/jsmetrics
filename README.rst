==================
jsmetrics: Jet-stream metrics and algorithms
==================
[Links]
----

This is ...

Inspired by xclim ()

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

DISCLAIMER
------------------
I have tried to replicate the various metrics based on the equations and details in the methodology as accurately as possible.

In some cases I have used a different dataset to the one used. 

This project is very much a work in progress, so contributors are very welcome. You  

Where you can find my working-out:
- I have included all of my working out in jupyter-notebooks available at: ... (warning: these notebooks have not been formatted nicely) 
- I am currently creating a verification notebook available at: ... where 


Table of the metrics
-------------------
See [jsmetrics/details_for_all_metrics.py] for specifications of each 
For their progress see [Project 1]

[TABLE HERE]


.. Documentation
.. -------------
.. The official documentation is at https://xclim.readthedocs.io/

.. Contributing
.. ------------
.. xclim is in active development and it's being used in production by climate services specialists.

.. * If you're interested in participating in the development of xclim by suggesting new features, new indices or report bugs, please leave us a message on the `issue tracker`_. There is also a chat room on gitter (|gitter|).

.. * If you would like to contribute code or documentation (which is greatly appreciated!), check out the `Contributing Guidelines`_ before you begin!

.. .. _issue tracker: https://github.com/Ouranosinc/xclim/issues
.. .. _Contributing Guidelines: https://github.com/Ouranosinc/xclim/blob/master/.github/CONTRIBUTING.rst


.. How to cite this library
.. ------------------------
.. If you wish to cite `xclim` in a research publication, we kindly ask that you use the bibliographical reference information available through `Zenodo`


Project To-Do's
-------------------
- ADD: cf_xarray (see: https://cf-xarray.readthedocs.io/en/latest/index.html)
- ADD: pint (see: https://pint.readthedocs.io/en/stable/)
- ADD: DOI
- LOOK INTO timing/benchmarking the metrics (maybe in seperate github repo)
- TO SOLVE: dealing with data from different sources (some sort of data translator module or maybe included in tests)
  - for example what if 'v' or 'v-wind' is passed to func instead of 'va' (answer: cf-xarray)
  - for example what if 'mbar' or 'model levels' instead of 'plev' (answer: pint)
- TO SOLVE: subsetting longitude if it wraps around 0-360
- CHECK: that methods using just U account for negative values (maybe need for abs() or not)