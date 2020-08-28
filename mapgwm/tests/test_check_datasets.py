""" Test the check_input_datasets routines. 
    As new routines are written, add to the import statement
    and then call appropriately
"""

import os
from mapgwm.check_input_datasets import check_headobs, compare_lists, get_header_length


def test_check_input_datasets(test_data_path):
    # input files for headobs
    data_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_stats_test.txt')
    metadata_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_meta_test.txt')
        
    check = check_headobs(data_file, metadata_file, os.path.join(os.path.join(os.getcwd(), 'mapgwm')))
    assert(check)
    

    
