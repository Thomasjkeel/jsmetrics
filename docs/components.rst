============================
Sub-components API reference
============================

What are the sub-components?
----------------------------
All statistics and algorithms in this package are built ontop of functions which we refer to as sub-components. 
These sub-component functions should have one role (e.g. to calculate atmospheric mass at a given kPa level), and should be easy to build from.
We have tried to make each sub-component as simple and single-use as possible, such that they could be used again to form future metrics implemented into the package. 
They are listed in full here.   

Jet statistic sub-components
----------------------------

.. automodule:: jsmetrics.metrics.jet_statistics_components
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Jet core algorithm sub-components
---------------------------------

.. automodule:: jsmetrics.metrics.jet_core_algorithms_components
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :noindex:

