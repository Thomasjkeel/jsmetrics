==============
Why jsmetrics?
==============



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
 

What are jet-streams?
---------------------
Jet streams are features of the atmospheric circulation that manifest as fast-flowing ribbons of air, usually around
8-12 km above the surface. 
They are generated and maintained in regions with extreme temperature gradients. These extreme gradients are produced
on the Earth by two major processes: (1) disturbances in the zonal mean-flow (known as eddy-driven processes) and (2)
conservation of angular momentum at the poleward edge of the Hadley Cell (known as thermally-driven processes).
In general, these processes create two major types of jets at a climatological scale:
1. the Polar Front Jet (PFJ) -- a deep and primarily eddy-driven feature
2. the Subtropical Jet (STJ) -- a shallow and primarily thermally-driven feature

.. image:: https://github.com/Thomasjkeel/jsmetrics/blob/write-docs/docs/_static/simple_jet_globe_diagram.jpeg
   :align: center
   :alt: Earth's two major jet streams

While it is often hard to seperate the two jets, the general strategy is to select 

They is no exact definition for what is and not a 'jet-stream' at any scale, as such differences in their measurement. 



.. Built from sub-components
.. ----------------------------
.. All statistics and algorithms in this package are built ontop of various one-purpose functions which we refer to as 'sub-components'. 
.. These sub-component functions should have one role (e.g. to calculate atmospheric mass at a given atmospheric level), and should allow yet to be added metrics an easier implementation.
 

