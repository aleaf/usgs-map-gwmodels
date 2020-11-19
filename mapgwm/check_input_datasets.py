"""Functions to test input files against those
   used for module testing by Travis CI.

   The automatic testing of the scripts uses standard
   input files in tests/data/ to make sure routines will
   run after changes are made.  These tests, however, do
   not check if there are changes to input file formats
   that can cause errors in the processing.  Functions in
   this package just test user specified input files against
   those used in the tests directory.  

   Authors should add a method for any datasets
   added to tests/data/
"""
from pathlib import Path
import numpy as np
import pandas as pd
import os
from mapgwm.headobs import preprocess_headobs, get_data


def get_header_length(sitefile, col0='SITE_BADGE'):
    with open(sitefile) as src:
        for i, line in enumerate(src):
            if '#' not in str(line) and col0 in str(line):
                return i

def compare_lists(list1, name1, list2, name2):
    """compare two lists and check for different or missing entries
       print out missing or different entries in list 2 compared to list 1

        Parameters
        ----------
        list1: list 
        name1: str - name to be printed in comparison
        list2: second list
        name2: str - name to be printed in comparison

        Returns
        -------
        passed: boolean,  returns True if files match
    """
    no_match = [x for x in list2 if x not in list1]
    missing = [x for x in list1 if x not in list2]
    passed = True
    if len(no_match) > 0:
        for x in no_match:
            print('{0} is in {2} but not in {1}'.format(x, name1, name2))
            passed = False
    elif len(missing) > 0:
        for x in missing:
            print('{0} is missing in {2} compared to {1}'.format(x, name1, name2))
            passed = False
    else:
        print('lists match: {0}, {1}'.format(name1, name2))
    return passed


def check_headobs_header(new_data_file, new_meta_file, path_to_tests='.'):
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
       (pass1 and pass2):  Boolean, False if files do not match
    """
       
    path = os.path.join(path_to_tests, 'tests', 'data', 'headobs')
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

    pass1 = compare_lists(testdata.columns.tolist(), 'Test Dataset', 
                  newdata.columns.tolist(), 'New Dataset')

    pass2 = compare_lists(testmeta.columns.tolist(), 'Test Metafile', 
                  newmeta.columns.tolist(), 'New Metafile')
    
    return (pass1 and pass2)


def check_preprocess_headobs_input(data_file, metadata_file, output_path='.'):
    test_data_path = Path(__file__).parent / 'tests/data'
    geographic_groups = [test_data_path / 'extents/CompositeHydrographArea.shp',
                         test_data_path / 'extents/MAP_generalized_regions.shp'
                         ]
    output_path = Path(output_path)
    outfile = output_path / 'preprocessed_monthly_output.csv'
    data, metadata = get_data(data_file=data_file, metadata_files=metadata_file)
    preproc_data, preproc_metadata = preprocess_headobs(data, metadata,
                                        head_data_columns=['head', 'last_head'],
                                        data_length_units='feet',
                                        active_area=test_data_path / 'extents/ms_delta.shp',
                                        source_crs=4269, dest_crs=5070,
                                        start_date='1998-04-01',
                                        geographic_groups=geographic_groups,
                                        geographic_groups_col='obsgroup',
                                        outfile=outfile)
    assert outfile.exists()
    assert Path(outfile.parent, outfile.stem + '_info.csv').exists()
    assert Path(outfile.parent, outfile.stem + '_info.shp').exists()
    assert np.all(preproc_data.columns ==
                  ['site_no', 'datetime', 'head', 'last_head', 'head_std', 'n', 'obsprefix'])
    assert not any(set(preproc_data.obsprefix).difference(preproc_metadata.obsprefix))
    assert not any({'site_no', 'x', 'y', 'screen_botm', 'screen_top',
                    'category', 'group'}.difference(preproc_metadata.columns))
    assert preproc_metadata['n'].dtype == np.integer
    # unit conversion was applied evenly
    preproc_data['head_diffs'] = np.abs(preproc_data['head'].values - preproc_data.last_head.values)
    preproc_data.sort_values(by='head_diffs', ascending=False, inplace=True)
    print(preproc_data.head(10))
    assert np.allclose(preproc_data['head'].values, preproc_data.last_head.values, atol=22)
    preproc_metadata['head_diffs'] = np.abs(preproc_metadata['head'].values - preproc_metadata.last_head.values)
    preproc_metadata.sort_values(by='head_diffs', ascending=False, inplace=True)
    print(preproc_metadata.head(10))
    assert np.allclose(preproc_metadata['head'].values, preproc_metadata.last_head.values, atol=22)

    # no negative open intervals
    assert not np.any((preproc_metadata.screen_top - preproc_metadata.screen_botm) < 0)


if __name__ == '__main__':
    data_file = os.path.join('tests', 'data', 'headobs', 'GW_monthly_stats_test.txt')
    metadata_file = os.path.join('tests', 'data', 'headobs', 'GW_monthly_meta_test.txt')

    check = check_headobs_header(data_file, metadata_file)
    if check:
        print('PASSED')
    else:
        print('FAILED')
