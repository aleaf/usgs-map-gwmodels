import os
import numpy as np
from mapgwm.headobs import preprocess_headobs, get_data


def test_preprocess_headobs(test_output_folder, test_data_path):
    # input files
    #data_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_stats_test.txt')
    data_file = test_data_path / 'headobs/GW_monthly_stats1990-01-01_2019-12-31.txt'
    #metadata_file = os.path.join(test_data_path, 'headobs', 'GW_monthly_meta_test.txt')
    metadata_file = test_data_path / 'headobs/GW_monthly_meta1990-01-01_2019-12-31.txt'

    # output
    outputfile = os.path.join(test_output_folder, 'preprocessed_monthly_output.csv')

    start_date = '1998-04-01'
    
    # areas of interest within model to break out as separate observation groups
    geographic_groups = [test_data_path / 'extents/CompositeHydrographArea.shp',
                         test_data_path / 'extents/MAP_generalized_regions.shp'
                         ]

    # read the data
    data_orig, metadata_orig = get_data(data_file, metadata_file)
    data, metadata = preprocess_headobs(data_orig, metadata_orig,
                                        head_data_columns=['head', 'last_head'],
                                        data_length_units='feet',
                                        active_area=os.path.join(test_data_path, 'extents/ms_delta.shp'),
                                        source_crs=4269, dest_crs=5070,
                                        start_date=start_date,
                                        geographic_groups=geographic_groups,
                                        geographic_groups_col='obsgroup',
                                        outfile=outputfile)
    assert os.path.exists(os.path.join(test_output_folder, 'preprocessed_monthly_output_info.shp'))
    assert os.path.exists(os.path.join(test_output_folder, 'preprocessed_monthly_output_info.csv'))
    assert np.all(data.columns ==
                  ['site_no', 'datetime', 'head', 'last_head', 'head_std', 'n', 'obsprefix'])
    assert not any(set(data.obsprefix).difference(metadata.obsprefix))
    assert not any({'site_no', 'x', 'y', 'screen_botm', 'screen_top',
                    'category', 'group'}.difference(metadata.columns))
    assert metadata['n'].dtype == np.integer
    # unit conversion was applied evenly
    assert np.allclose(data['head'].values, data.last_head.values, rtol=0.1)
    assert np.allclose(metadata['head'].values, metadata.last_head.values, rtol=0.1)
    if data_orig.head_std.any() and data.head_std.any():
        assert np.allclose(np.nanmean(data_orig.head_std)/np.nanmean(data.head_std), 3.28, rtol=0.1)

    # no negative open intervals
    assert not np.any((metadata.screen_top - metadata.screen_botm) < 0)
