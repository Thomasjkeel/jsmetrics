=======
History
=======

0.2.9 (2025-01-16)
-------------------------
* Fix ongoing numpy dependency issue and update requirement for numpy version so it is less strict
* Remove warnings.RankWarning fix for quadratic function. This should fix Cattiaux

0.2.8 (2024-12-17)
-------------------------
* Fix scipy.optimize.least_sq depreciation warning
* Fix xarray groupby squeeze depreciation warning
* Fix numpy depreciation warning

0.2.7 (2024-09-24)
-------------------------
* Fix scipy.ndimage maximum filter error
* Corrects Kuang et al. 2014 and Schiemann et al. 2009 methods

0.2.6 (2024-09-23)
-------------------------
* Add expand_dims and squeeze to time dimensions of 'run_jet_core_and_region_algorithm_on_one_day' to fix manney_et_al_2011


0.2.4 (2024-09-23)
-------------------------
* Fix scipy.ndimage depreciation warning
* Change argmax to nanargmax
* update max numpy version to 1.26 so 'np.RankWarning' still exists (in future prepare for np v2.0)
* Make sure requirement files have same versions


0.2.3 (2024-05-10)
-------------------------
* Add example of assigning units to plev in for docstrings in data_utils
* Add step-security-bot and dependenabot for OSSF scorecard
* Update dependencies and workflows


0.2.3-alpha (2024-04-08)
-------------------------
* Add example of assigning units to plev in data_utils


0.2.2 (2024-03-15)
-------------------------
* Overhaul Kuang et al. 2014 method to work with plev-lat center check rather than lat-long jet center check.
* Use maximum filter from scipy ndimage for diagonal check from 'find_local_maxima_in_2d_dataarray_with_diagonals' used in M11 and K14
* Also makes a change to Manney et al. 2011 method if using diagonal check (should not affect results as it was working before)


0.2.2-beta (2024-03-13)
-------------------------
* Fix error when checking for plev units in data_utils
* Reduce README content


0.2.1 (2023-12-01)
-------------------------
* Update Manney et al. 2011 method to include 'check_diagonals' parameter
* Update Manney et al. 2011 method to select maximum windspeed if no cores found above threshold
* Add 'slice_array_by_index_breaks' to data utils
* Update examples of use and quick start table in docs
* Add check for plev units to data_utils
* Update Archer & Calidera to more properly check for whether data plev coords are Pa or hPa.
* Fix typos and clarify writing in docs


0.2.0 (2023-10-12)
-------------------------
* Fix CI and Github actions
* Finish overhaul the ReadTheDocs documentation in the package
* Confirm completed and to verify metrics
* Update Schiemann et al. 2009 method to run faster and introduce u_threshold parameter default
* Add more complete inline 'Steps' to each metric in the package for readability.
* Remove Screen & Simmonds and Local Wave Activity leftover code
* Update Manney et al. 2011 further with padding to capture cores near the edge of dataset

0.1.7 (2023-10-10)
-------------------------
* Changes:
* **MAJOR:** Rename to Barnes & Polvani 2014 and update tests and all mentions of
* Update Cattiaux method to deal with Matplotlib depreciation warning about using 'allsegs' or 'collections' in contour plot
* Add tests for new spatial utils method 'seperate_one_contour_into_line_segments'
* Add readthedocs.yml

0.1.7-beta (2023-10-05)
-------------------------
* Rename to Barnes & Polvani 2014 and update tests and all mentions of

0.1.7-alpha (2023-09-19)
-------------------------
* Update Cattiaux method to deal with Matplotlib depreciation warning about using 'allsegs' or 'collections' in contour plot
* Add tests for new spatial utils method 'seperate_one_contour_into_line_segments'

0.1.6 (2023-09-14)
-------------------------
* Overhaul the ReadTheDocs documentation in the package
    * Write a more detailed description of each metric
    * Add examples and implementation notes to each metric
    * Update listing in 'details_for_all_metrics.py'
    * Update the 'Why jsmetrics' section with notes about what jets are and quick start help
    * Add some Examples of Use for each type of metric
* Add ws_threshold parameter to Pena-Ortiz method
* Rename Pena-Ortiz method to 'pena_ortiz_et_al_2016' (as earlier paper)

0.1.6-alpha (2023-08-19)
-------------------------
* Fix Manney et al. 2011 implementation
* Correct Manney 2011 method
* Move old method to new metric: 'jet_core_identification_algorithm'
* Update docs for Manney 2011 and sub-components

0.1.5-beta (2023-08-16)
-------------------------
* Update Schiemann method with new variable name (jet occurence), docstring and changes to sub-component function names
* Add ws_threshold parameter to Schiemann
* Update methods that work on one time unit of data and add squeeze method to properly deal if time not in dims
* Add basic outline of jet core algorithm docstrings

0.1.5-alpha (2023-08-15)
-------------------------
* Begin overhaul of ReadTheDocs documentation
* Add notes and example to Koch et al. 2006 metric
* Rename variable returned by K06 to jet_events_ws

0.1.4 (2023-08-02)
-------------------------
* Add new metric to package: Zappa et al. 2018 (This method builds on Ceppi et al. 2018)
* Raise KeyError if no time coordinate is passed to a given metrics

0.1.4-alpha (2023-07-21)
-------------------------
* add KeyError raise if no time coordinate is passed to various metrics

0.1.3 (2023-07-07)
-------------------------
* Add "method='nearest'" to jet statistics and core algorithms for cases when coords cannot be represented within float precision range.

0.1.2 (2023-06-06)
-------------------------
* Fix  Barnes & Polvani 2013 to better deal when min max jet lat is at edge data
* Add check for NoLeapDatetime

0.1.2-alpha (2023-05-27)
-------------------------
* Add check for NoLeapDatetime

0.1.1 (2023-05-26)
-------------------------
* Fix Woollings et al. 2010 and filter windows to use day timeunits for window to stop it removing too much data.
* Add data util function to add number of days to 360Day Datetime type

0.1.1-beta (2023-04-07)
-------------------------
* add parameter for Kerr et al. 2020
* Add Ceppi et al jet speed adaptation from Screen et al. 2022
* Add fix for sort_xarray_data_coords so it works when only one coord value in coordinate (i.e. so each metric can work when only one longitude)
* Supress warning for quadratic func


0.1.1-alpha (2023-03-31)
-------------------------
* Add fix for Kuang to run when there is no time dim
* Add fix for BP15 to except errors where all nan data
* Add warning for BS17 when more than 10 days resolution

0.1.0 (2023-01-22)
-------------------------
* MAJOR UPDATE: re-organise the structure of the package into core, metrics and utils
* rename jet statistics, waviness metrics and jet core algorithm files
* add wrappers to check data is xarray and is sorted in descending order (in core/check_data.py)
* move waviness metrics to new file
* Update appropriate tests

0.0.19-alpha (2022-12-21)
-------------------------
* Update JetStreamOccurenceAndCentreAlgorithm to skip longitude values outside lon range in data
* Make changes to work with Shapely version 1.8/2.0. Means changes to Cattiaux et al. 2016

0.0.18 (2022-11-23)
-------------------------
* update fitted parabola func for Barnes & Polvani 2015
* Add Blackport & Fyfe 2022
* update Barnes & Simpson 2017 to drop all NaN slices
* update to check for more than one time step for time groupby methods
* add test to check all metrics when input is one time step

0.0.17 (2022-11-13)
-------------------------
* add try and except for Grise & Polvani 2017 to account for missing vals

0.0.16 (2022-11-09)
-------------------------
* skipna=True for calc_latitude_and_speed_where_max_ws
* Barnes and Simpson mean over longitude for jet lat

0.0.15 (2022-11-09)
-------------------------
* rename max_lat_0.01 to jet_lat for Grise & Polvani 2017
* Fix get_3_latitudes_and_speed_around_max_ws to work with isel around lon
* Fix barnes polvani parabola to deal with nan values

0.0.14 (2022-11-09)
-------------------------
* add plev mean to Bracegirdle

0.0.14-alpha (2022-10-25)
-------------------------
* update Pena Ortiz so that it returns monthyear and by day local wind maxima
* remove make_empty_local_wind_maxima_data func
* Fix CI
* Add millibars to get_all_hPa_list

0.0.13 (2022-10-19)
-------------------------
* fix workflow for publish to PyPi and TestPyPi

0.0.12 (2022-10-19)
-------------------------
* fix kuang to work for southern hemisphere as well
* add workflow for publish to PyPi

0.0.12-alpha (2022-10-18)
-------------------------
* Update calc_latitude_and_speed_where_max_ws to use numpy methods
* Fix Barnes and Simpson 2017 method so it runs on each longitude

0.0.11 (2022-09-15)
-------------------------
* Update and fix the JetStreamOccurenceAndCentreAlgorithm method for Kuang
* Change LICENSE
* Upload to Zenodo

0.0.10 (2022-08-21)
-------------------
* First release to pypi
* Clean up rst docs

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
