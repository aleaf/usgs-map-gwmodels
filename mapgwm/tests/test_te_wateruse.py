import os
import numpy as np
import pytest
from mapgwm.te_wateruse import preprocess_te_wateruse, read_te_water_use_spreadsheet


@pytest.mark.parametrize(('xlsx_file,sheet_name,skiprows,date,'
                          'x_coord_col,y_coord_col,q_col,site_no_col,'
                          'site_name_col,source_name_col,source_code_col'),
                         (('swuds/2010_Thermo_Model_Estimates.xlsx', 'Report_table_UPDATED', 2, '2010',
                           'Longitude, decimal degrees',
                           'Latitude, decimal degrees',
                           'USGS-Estimated Annual Withdrawal (Mgal/d)',
                           'PLANT CODE', 'PLANT NAME', 'NAME OF WATER SOURCE', 'USGS WATER SOURCE CODE'
                           ),
                          ('swuds/2015_TE_Model_Estimates_lat.long_COMIDs.xlsx', '2015_ANNUAL_WD_CU', 0, '2015',
                            'LATITUDE', 'LONGITUDE', 'WITHDRAWAL', 'EIA_PLANT_ID', 'PLANT_NAME',
                           'NAME_OF_WATER_SOURCE', 'WATER_SOURCE_CODE'
                           )
                          ))
def test_read_te_water_use_spreadsheet(test_data_path, xlsx_file, sheet_name, skiprows,
                                       date, x_coord_col, y_coord_col, q_col, site_no_col,
                                       site_name_col, source_name_col, source_code_col):
    xlsx_file = test_data_path / xlsx_file
    df = read_te_water_use_spreadsheet(xlsx_file, date='2010', x_coord_col=x_coord_col,
                                       y_coord_col=y_coord_col,
                                       q_col=q_col, site_no_col=site_no_col,
                                       site_name_col=site_name_col, source_name_col=source_name_col,
                                       source_code_col=source_code_col,
                                       sheet_name=sheet_name, skiprows=skiprows)
    assert np.all(df.columns ==
                  ['site_no', 'start_datetime', 'x', 'y', 'q', 'site_name', 'source_name'])
    assert 'datetime64' in df.start_datetime.dtype.name


@pytest.fixture(scope='function')
def te_data(test_data_path, test_output_folder):
    df2010 = read_te_water_use_spreadsheet(test_data_path / 'swuds/2010_Thermo_Model_Estimates.xlsx',
                                           date='2010',
                                           x_coord_col='Longitude, decimal degrees',
                                           y_coord_col='Latitude, decimal degrees',
                                           q_col='USGS-Estimated Annual Withdrawal (Mgal/d)',
                                           site_no_col='PLANT CODE',
                                           site_name_col='PLANT NAME',
                                           source_name_col='NAME OF WATER SOURCE',
                                           source_code_col='USGS WATER SOURCE CODE',
                                           sheet_name='Report_table_UPDATED', skiprows=2)

    df2015 = read_te_water_use_spreadsheet(
        test_data_path / 'swuds/2015_TE_Model_Estimates_lat.long_COMIDs.xlsx',
        date='2015',
        x_coord_col='LONGITUDE',
        y_coord_col='LATITUDE',
        q_col='WITHDRAWAL',
        site_no_col='EIA_PLANT_ID',
        site_name_col='PLANT_NAME',
        source_name_col='NAME_OF_WATER_SOURCE',
        source_code_col='WATER_SOURCE_CODE',
        sheet_name='2015_ANNUAL_WD_CU', skiprows=0)
    df = df2010.append(df2015)
    return df


@pytest.fixture(scope='function')
def preprocessed_te_data(te_data, test_data_path, test_output_folder):
    df = te_data
    outfile = test_output_folder / 'preprocessed_te_data.csv'

    estimated_production_zone_top = test_data_path / 'swuds/rasters/pz_MCAQP_top.tif'
    estimated_production_zone_botm = test_data_path / 'swuds/rasters/pz_MCAQP_bot.tif'

    results = preprocess_te_wateruse(df,  # dem=dem, dem_units='feet',
                                     start_date='2008-01-01',
                                     end_date='2017-12-31',
                                     active_area=os.path.join(test_data_path, 'extents/ms_delta.shp'),
                                     estimated_production_zone_top=estimated_production_zone_top,
                                     estimated_production_zone_botm=estimated_production_zone_botm,
                                     estimated_production_surface_units='feet',
                                     source_crs=4269,
                                     dest_crs=5070,
                                     data_volume_units='mgal',
                                     model_length_units='meters',
                                     outfile=outfile)
    return results


def test_preprocess_te_wateruse(preprocessed_te_data, test_data_path, test_output_folder):

    outfile = test_output_folder / 'preprocessed_te_data.csv'
    results = preprocessed_te_data
    assert np.all(results.screen_top > results.screen_botm)
    assert np.all(results.columns[:8] ==
                  ['site_no', 'start_datetime', 'x', 'y', 'screen_top', 'screen_botm', 'q', 'geometry'])
    assert outfile.exists()
    assert outfile.with_suffix('.shp').exists()
    check_cols = ['q', 'screen_top', 'screen_botm', 'x', 'y', 'start_datetime', 'site_no']
    # discharges must be negative!
    assert results['q'].sum() < 0
    for col in check_cols:
        assert not results[col].isnull().any()