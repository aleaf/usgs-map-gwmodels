import os
import numpy as np
from mapgwm.headobs import preprocess_headobs, get_data


def test_preprocess_headobs(test_output_folder, test_data_path):
    # input files
    data_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_stats_test.txt')
    metadata_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_meta_test.txt')

    # output
    outputfile = os.path.join(test_output_folder, 'preprocessed_monthly_output.csv')

    start_date = '1998-04-01'
    
    # areas of interest within model to break out as separate observation groups
    aoi = {'DeltaAOI': os.path.join(test_data_path, 'extents', 'CompositeHydrographArea.shp')
           }

    # read the data
    data, metadata = get_data(data_file, metadata_file)
    data, metadata = preprocess_headobs(data, metadata,
                                        data_length_units='feet',
                                        active_area=os.path.join(test_data_path, 'extents/MERAS_Extent.shp'),
                                        src_crs=4269, dest_crs=5070,
                                        aoi=aoi, start_date=start_date,
                                        outfile=outputfile)
    assert os.path.exists(os.path.join(test_output_folder, 'preprocessed_monthly_output_info.shp'))
    assert os.path.exists(os.path.join(test_output_folder, 'preprocessed_monthly_output_info.csv'))
    assert np.all(data.columns ==
                  ['site_no', 'datetime', 'head_m', 'last_head_m', 'head_std_m', 'n', 'obsprefix'])
    assert metadata.index.name == 'site_no'
    assert not any(set(data.obsprefix).difference(metadata.obsprefix))
    assert not any({'x_5070', 'y_5070', 'screen_botm', 'screen_top',
                    'category', 'group'}.difference(metadata.columns))
