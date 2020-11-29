import os
import numpy as np
import pytest
from mfsetup.units import convert_length_units, convert_volume_units
from mapgwm.swuds import Swuds, preprocess_swuds


@pytest.fixture(scope='function')
def preprocessed_swuds_data(test_data_path, test_output_folder):

    swuds_input = os.path.join(test_data_path, 'swuds', 'withdrawals_test_dataset.xlsx')
    worksheet = 'LMG-withdrawals-2000-2018'
    outfile = test_output_folder / 'preprocessed_swuds_data.csv'

    dem = os.path.join(test_data_path, 'rasters', 'dem_mean_elevs.tif')

    production_zones = {'mrva': (test_data_path / 'swuds/rasters/pz_MCAQP_top.tif',
                                 test_data_path / 'swuds/rasters/pz_MCAQP_bot.tif',
                                 ),
                        'mcaq': (test_data_path / 'swuds/rasters/pz_MCAQP_top.tif',
                                 test_data_path / 'swuds/rasters/pz_MCAQP_bot.tif',
                                 'feet'),
                        'lcaq': (test_data_path / 'swuds/rasters/pz_LCAQP_top.tif',
                                 test_data_path / 'swuds/rasters/pz_LCAQP_bot.tif')
                        }

    results = preprocess_swuds(swuds_input, worksheet,
                               dem=dem, dem_units='feet',
                               start_date=None,
                               end_date=None,
                               active_area=os.path.join(test_data_path, 'extents/ms_delta.shp'),
                               production_zones=production_zones,
                               source_crs=4269,
                               dest_crs=5070,
                               data_length_units='feet',
                               data_volume_units='mgal',
                               model_length_units='meters',
                               outfile=outfile)
    return results


def test_preprocess_swuds(preprocessed_swuds_data, test_data_path, test_output_folder):

    outfile = test_output_folder / 'preprocessed_swuds_data.csv'
    results = preprocessed_swuds_data
    assert results.data_length_units == 'feet'
    assert np.allclose(convert_length_units(results.data_length_units,
                                            results.model_length_units), 0.3048)
    assert np.allclose(convert_volume_units(results.data_volume_units,
                                            results.model_length_units), 3785.4, atol=0.02)
    assert outfile.exists()
    assert outfile.with_suffix('.shp').exists()
    check_cols = ['q', 'screen_top', 'screen_bot', 'x', 'y', 'start_datetime', 'site_no']
    # discharges must be negative!
    assert results.df['q'].sum() < 0
    for col in check_cols:
        assert not results.df[col].isnull().any()
