=============
API reference
=============

Metric sub-components
#####################

What are the sub-components?
----------------------------
All statistics and algorithms in this package are built ontop of functions which we refer to as sub-components.
These sub-components (should) have one role (e.g. to calculate atmospheric mass at a given kPa level), and can be used to build future metrics implemented into the package upon.
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


Utility functions
#################
All functions used in package that provide generic utility (such as for xarray compatibility) are listed here, and are seperate from metric sub-components.

Data Utils
------------------------

.. automodule:: jsmetrics.utils.data_utils
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Spatial Utils
------------------------

.. automodule:: jsmetrics.utils.spatial_utils
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Windspeed Utils
------------------------

.. automodule:: jsmetrics.utils.windspeed_utils
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :noindex:
