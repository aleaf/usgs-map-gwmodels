import os
from mapgwm.headobs import preprocess_headobs, get_data, get_active_area


def test_preprocess_headobs(test_output_folder, test_data_path):
    # input files
    data_file = os.path.join(test_data_path, 'headobs/GW_monthly_stats1990-01-01_2019-12-31.txt')
    metadata_file = os.path.join(test_data_path, 'headobs/GW_monthly_meta1990-01-01_2019-12-31.txt')

    # output
    outputfile = os.path.join(test_output_folder, 'preprocessed_monthly_output.csv')

    start_date = '1998-04-01'
    
    # areas of interest within model to break out as separate observation groups
    aoi = {'DeltaAOI': os.path.join(test_data_path, 'extents/CompositeHydrographArea.shp')
           }
    
    # area of observations to process (discard observations outside of this area)
    active_area = get_active_area(os.path.join(test_data_path, 'extents/MERAS_Extent.shp'),
                                  name_col='FORMATION',
                                  buffer=10000.)

    # read the data
    data, metadata = get_data(data_file, metadata_file)
    preprocess_headobs(data, metadata, data_length_units='feet',
                       active_area=active_area,
                       aoi=aoi, start_date=start_date,
                       outfile=outputfile)
    assert os.path.exists(os.path.join(test_output_folder, 'preprocessed_monthly_output_info.shp'))
    assert os.path.exists(os.path.join(test_output_folder, 'preprocessed_monthly_output_info.csv'))
