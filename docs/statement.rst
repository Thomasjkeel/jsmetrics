==============
Why jsmetrics?
==============
The planet see
 


Quick start
-----------
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
 

What are jet streams?
---------------------
Jet streams are features of the atmospheric circulation that manifest as fast-flowing ribbons of air, usually around
8-12 km above the surface. 
They are generated and maintained in regions with extreme temperature gradients. These extreme gradients are produced
on the Earth by two major processes: (1) disturbances in the zonal mean-flow (known as eddy-driven processes) and (2)
conservation of angular momentum at the poleward edge of the Hadley Cell (known as thermally-driven processes).

In general, these processes create two major types of jets at a climatological scale in each Hemisphere (see Figure 1):

   1. the Polar Front Jet (PFJ) -- a deep and primarily eddy-driven feature
   2. the Subtropical Jet (STJ) -- a shallow and primarily thermally-driven feature

.. figure:: _static/images/simple_jet_globe_diagram.jpeg
   :align: center
   :alt: Earth's two major jet streams

   Figure 1. Idealised view of the planet's jet streams

Figure 1 shows a idealised version of the jet streams -- clearly seperated and flowing circumglobal west-to-east fashion.
As you can imagine, in reality, the location, strength and direction of a given jet stream is not well defined at any scale.
They also exhibit fairly strong seasonality (generally moving closer to the Equator in colder months, and closer to the poles
in warmer ones). To see what we mean, we reccomend having a play with `Earth null school <https://earth.nullschool.net/#2021/02/15/1700Z/wind/isobaric/250hPa/orthographic=-91.82,32.12,310>`_
for one view of how jet stream-like features manifest on the planet (i.e. at 250-850 hPa).

The complexity in their structure, and lack of strong definition (they are essential just 'atmospheric phenomena') means that a
vast range of metrics, statistics and algorithms have been employed to identify and characterise different aspects of them in
atmospheric data. With *jsmetrics*, we have tried to include as many of the most common methods used to characterise jet streams
as possible in the hope that this would help researchers reconcile information about them and allow for a more quantitative
comparison of their differences and impact on trends and changes shown to the jet streams.

*I am still writing this section, so please email me if you have some suggestions or feedback.*


.. Built from sub-components
.. ----------------------------
.. All statistics and algorithms in this package are built ontop of various one-purpose functions which we refer to as 'sub-components'. 
.. These sub-component functions should have one role (e.g. to calculate atmospheric mass at a given atmospheric level), and should allow yet to be added metrics an easier implementation.
 

