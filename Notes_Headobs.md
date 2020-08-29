Notes on headobs.py
===================

Notes on collecting headobs scripts from meras6, cache, grandprairie, delta 
and central sands repositories and making a common *headobs.py* script that 
can be used in the preprocessing for any of the MAP models.  The script 
reads in groundwater head observation data as processed by Will Asquith's 
routines.

Testing
-------

The script *test_headobs.py* in the *tests* subdirectory should produce
a csv and shapefile of monthly groundwater heads from the input file specified
in the test.

Original
--------

There was a headobs.py script in the directory already.  It was compared 
to the other scripts and updated as described below.


Cache
-----

- Preprocess_headobs.py for the Cache model has a function *assign_wells_to_layers()*, 
  but it depends on having a modflow model object.  The x,y location of the well is used to 
  assign row and column; and then the appropriate model layer is identified given information
  about the well screen or open interval of the well.  This seems like a modflow-setup
  step, so it was not kept in the *headobs.py* script.  
  
- Imported zipfile and datetime as dt -- neither seemed to be needed

- In *get_active_area()*, the original version had df.index = df[name_col] commented
  out.  But, with that commenting, the name_col parameter is not used.  I uncommented
  in the *headobs.py* script to match Cache.


Grand Prairie
-------------

The scripts folder has a notebook *observations.ipynb* that uses
Leaf's *pydrograph* package to pull data from NWIS and processes
that information.  


Delta
-----

Does not have a separate headobs preprocessing script.


Meras3
------

This script also imported zipfile and datetime.  It uses an older workflow
where the input datasets were read from csv files that had been zipped.  
It also has the *assign_wells_to_layers()* function which is not used
for preprocessing.  The script in *mapgwm* was retained as the latest version.



  
  


