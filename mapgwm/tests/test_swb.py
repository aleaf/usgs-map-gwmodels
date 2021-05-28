import os
from mapgwm.swb import get_monthly_values
from mfsetup import load_modelgrid


def test_preprocess_swb(test_output_folder, test_data_path):
    # input files
    data_file = os.path.join(test_data_path, 'swb', 'cache_input_August2020.nc')

    # output
    outputfile = os.path.join(test_output_folder, 'swb_monthly_means.nc')

    # get area of SWB to process
    active_grid = load_modelgrid(os.path.join(test_data_path, 
                                             'extents', 
                                             'cache250_grid.json'))

    # get the monthly means
    get_monthly_values(data_file, outfile=outputfile,
                      filter=active_grid.bounds, check_results=True)
    assert os.path.exists(os.path.join(test_output_folder, 'swb_monthly_means.nc'))
    
