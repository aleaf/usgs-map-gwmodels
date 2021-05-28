""" Test the check_input_datasets routines. 
    As new routines are written, add to the import statement
    and then call appropriately
"""

import os
from mapgwm.check_input_datasets import check_headobs_header, compare_lists, get_header_length


def test_check_input_datasets(test_data_path):
    # input files for headobs
    #data_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_stats_test.txt')
    #metadata_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_meta_test.txt')
    data_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_stats1990-01-01_2019-12-31.txt')
    metadata_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_meta1990-01-01_2019-12-31.txt')
    assert os.path.exists(data_file), data_file
    
    # check_headobs_header returns True if both data_file and metadata_file match
    check = check_headobs_header(data_file, metadata_file, test_data_path)
    assert(check)
    

    
