import os
import numpy as np
import pandas as pd
import pytest
from mapgwm.iwum import preprocess_iwum_pumping, plot_iwum_output


@pytest.fixture
def ncfile(test_data_path):
    return test_data_path / 'iwum/irr_wu_1km_1999_to_2017.nc'


def test_preprocess_iwum_pumping(test_data_path, test_output_folder, ncfile):
    production_zone_top_raster = test_data_path / 'iwum/est_prod_zone_top.tif'
    production_zone_botm_raster = test_data_path / 'iwum/est_prod_zone_botm.tif'
    outfile = test_output_folder / 'iwum_m3.csv'
    results = preprocess_iwum_pumping(ncfile,
                                      start_date=None,
                                      end_date=None,
                                      active_area=os.path.join(test_data_path, 'extents/ms_delta.shp'),
                                      estimated_production_zone_top=production_zone_top_raster,
                                      estimated_production_zone_botm=production_zone_botm_raster,
                                      flux_variable='value',
                                      nc_crs=5070,
                                      nc_length_units='meters',
                                      estimated_production_surface_units='meters',
                                      model_length_units='meters',
                                      outfile=outfile)
    assert np.all(results.screen_top >= results.screen_botm)
    assert np.all(results.columns ==
                  ['site_no', 'x', 'y', 'screen_top', 'screen_botm', 'start_datetime',
                   'end_datetime', 'q'])
    assert outfile.exists()
    check_cols = ['q', 'x', 'y', 'start_datetime', 'end_datetime',
                  'screen_top', 'screen_botm', 'site_no']
    for col in check_cols:
        assert not results[col].isnull().any()


def test_plot_iwum_output(ncfile, test_output_folder):
    results = plot_iwum_output(ncfile,
                               flux_variable='value',
                               output_path=test_output_folder)
    ftime = pd.Timestamp(os.path.getmtime(ncfile), unit='s')
    outfile = test_output_folder / (ncfile.stem + '_{}.pdf'.format(ftime.strftime('%Y-%m-%d')))
    assert outfile.exists()
