"""Functions to test input files against those
   used for module testing by Travis CI.

   The automatic testing of the scripts uses standard
   input files in tests/data/ to make sure routines will
   run after changes are made.  These tests, however, do
   not check if there are changes to input file formats
   that can cause errors in the processing.  Functions in
   this package just test user specified input files against
   those used in the tests directory.  

   Authors should add a check_ method for any datasets
   added to tests/data/
"""


import pandas as pd
import os


def get_header_length(sitefile, col0='SITE_BADGE'):
    
    with open(sitefile) as src:
        for i, line in enumerate(src):
            if '#' not in str(line) and col0 in str(line):
                return i


def test_headobs(new_data_file, new_meta_file):
    """headobs.py reads in data from Will Asquith as two csv files
       datafile - csv with observed monthly heads
       metadata - csv with site characteristics

       test_headobs points to 
       mapgwm/tests/data/headobs/GW_monthly_stats_test.txt
       mapgwm/tests/data/headobs/GW_monthly_meta_test.txt

       Which have the top 1000 lines of datasets from August 2020.

       Parameters
       ----------
       new_data_file: path to data file to test
       new_meta_file: path to metadata file to test


       Returns
       -------
       None:  will throw assert error if data files do not match
    """
       
    path = os.path.join('tests', 'data', 'headobs')
    data_file = os.path.join(path, 'GW_monthly_stats_test.txt')
    metadata_file = os.path.join(path, 'GW_monthly_meta_test.txt')

    data_skiprows = get_header_length(data_file)
    testdata = pd.read_csv(data_file, sep='\t', skiprows=data_skiprows)
    data_skiprows = get_header_length(metadata_file)
    testmeta = pd.read_csv(metadata_file, sep='\t', skiprows=data_skiprows)

    data_skiprows = get_header_length(new_data_file)
    newdata = pd.read_csv(new_data_file, sep='\t', skiprows=data_skiprows)
    data_skiprows = get_header_length(new_meta_file)
    newmeta = pd.read_csv(new_meta_file, sep='\t', skiprows=data_skiprows)

    assert (testdata.columns.tolist().sort() == newdata.columns.tolist().sort())
    assert (testmeta.columns.tolist().sort() == newmeta.columns.tolist().sort())
    print('headobs data files matched')


if __name__ == '__main__':
    path = os.path.join('tests', 'data', 'headobs')
    data_file = os.path.join(path, 'GW_monthly_stats_test.txt')
    metadata_file = os.path.join(path, 'GW_monthly_meta_test.txt')
    test_headobs(data_file, metadata_file)