====================
Why use jsmetrics?
====================


What are the sub-components?
----------------------------
All statistics and algorithms in this package are built ontop of various one-purpose functions which we refer to as sub-components. 
These sub-component functions should have one role (e.g. to calculate atmospheric mass at a given atmospheric level), and should be easy to build from.
This, it is hoped, will make future additions of metrics easier to implement to this package, as well as allowing the users to tweak aspects of the current metrics in the package (e.g. allowing the user to swap out a given smoothing function)

Quick help
------------------------
.. table::
   :align: left
   :widths: auto
   
   ======================================================= ===============================================
   I would like to...                                      Reccomendation 
   ======================================================= ===============================================
   Know the average position of the jet stream             Use a jet latitude statistic
   Know the average speed of the jet stream                Use a jet latitude statistic
   Make a map of the jet stream                            Use a jet core algorithm to   
   ...                                                     ...
   ======================================================= ===============================================
 

