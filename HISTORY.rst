=======
History
=======


0.0.10 (2022-08-21)
-------------------
* First release to pypi
* Clean up rst docs
* remove dependency versions

0.0.9 (2022-08-16)
------------------
* Finish tests
* Remove TODOs
* Outline metric_verification notebooks
* Improve docs

0.0.8 (2022-07-18)
------------------
* Format the readme
* seperate metrics into metrics and algorithms
* Reorder and write better docstrings for the utils files 
* Update year on LICENSE 

0.0.7-beta (2022-06-30)
-----------------------
* swap 'plev' and 'lat' in manney_et_al_2011 method so that it groups cores better
* rename 'sinouisity' to 'sinuosity'

0.0.7-alpha (2022-06-10)
------------------------
* update spatial_utils with lazy method for guessing bounds and assuming a regular grid (func is "_standardise_diffs_by_making_all_most_common_diff")
* update Pena-Ortiz method to seperate into subtropical and polar front jet
* remove prints from windspeed utils
* rename bp13 jet lat 

0.0.6 (2022-06-09)
------------------
* add Barnes & Polvani 2015 
* add Kerr et al. 2020
* add nearest method function to general utils
* Speed up Ceppi and fix integration method within (still need to verify)
* Add spatial utils for grid cell m2 method

0.0.6-beta (2022-05-31)
-----------------------
* Fix 'get_latitude_and_speed_where_max_ws_at_reduced_resolution' with check for np.nans

0.0.6-alpha (2022-05-25)
------------------------
* add Barnes & Polvani 2013
* Fix 'get_latitude_and_speed_where_max_ws' so it can take one value 
* Fix Barnes & Simpson 2017 and Woollings et al. 2010 and change name of col
* Fix Barnes & Polvani neighbouring lats  and speed 

0.0.5 (2022-05-23)
------------------
* add Barnes & Simpson 2017 
* Update 'get_latitude_and_speed_where_max_ws' function 
* Update calc_mass_weighted wind 

BIG CHANGES
^^^^^^^^^^^
* Change the 'get_latitude_and_speed_where_max_ws' function to take abs() max -> will mean that negative u-wind values can be considered the jet lat


0.0.5-beta (2022-05-03)
-----------------------
* update Woollings et al. 2010 with seasonal cycle
* update metric details dict with 'plev_units' argument 
* fix archer and caldiera call to mass weighted ws (STILL TODO: better plev understanding)

0.0.5-alpha (2022-04-24)
------------------------
* add metric verification notebooks 

0.0.4-beta (2022-02-09)
-----------------------
* add description, name and DOI to metric details dict

0.0.4-alpha (2022-01-26)
------------------------
* remove Docker
* remove get data scripts

0.0.3-gamma (2022-01-14)
------------------------
* remove python 3.6 compatibility
* update environment yml (still broken)

0.0.3-beta (2022-01-14)
-----------------------
* Use real part from fourier filter to Woollings and its tests

0.0.3-alpha (2022-01-14)
------------------------
* Remove main and experiment related files (moved to another directory so this one is cleaner)

0.0.2 (2022-01-10)
------------------
* First release on github

0.0.2-beta (2022-01-10)
-----------------------

* Add docstrings to all metrics and sub-components

0.0.2-alpha (2022-01-04)
------------------------

* Add docstrings to Archer & Calidera metric

0.0.1 (2022-01-04)
------------------

* Allow jsmetric to call jetstream_metrics and utils

0.0.1-beta (2021-12-30)
-----------------------

* Add currently existing metrics
