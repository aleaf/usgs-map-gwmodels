Philosophy
==========

Motivation
----------
The motivation for this project is bring together common scripts used in 
the development of the regional and inset MAP models to a single package.
That way, if an input dataset from some other part of the MAP project is
changed, the necessary edits to the processing script can be made once.
If common scripts exist for each of the insets and the regional model, then
those edits would have to be made for each of them, increasing the chance
for error.


What usgs-map-gwmodels does
-------------------------------------------
Functions in this package can be used to help develop preprocess data for
MAP models.  The package can be used along with *modflow-setup* to develop
models from the various data sources.  Post-processing scripts could also
be put into this package, but they have not been added yet.

What usgs-map-gwmodels doesn't do
--------------------------------------------------------
The most important limitation is that this package is not designed to
provide set up functions that are specific to only one of the regional models.
Model-specific functions or scripts should be part of the client code (python
script or jupyter notebooks) that use *usgs-map-gwmodels* and *modflow-setup*
or user maintained packages. 