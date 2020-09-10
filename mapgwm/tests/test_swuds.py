import os
from mapgwm import swuds


def test_preprocess_swuds(test_output_folder, test_data_path):
   
    # output
    outcsv = os.path.join(test_output_folder, 'withdrawals_test.csv')
    outputshp = os.path.join(test_output_folder, 'swuds_select_wells.shp')
    processed_csv = os.path.join(test_output_folder, 'processed_swuds.csv')

    # required data
    swuds_input = os.path.join(test_data_path, 'swuds', 'withdrawals_test_dataset.xlsx')
    worksheet = 'LMG-withdrawals-2000-2018'

    # required rasters and shapefile
    dem = os.path.join(test_data_path, 'extents', 'dem_mean_elevs_1000.tif')
    mc_top = os.path.join(test_data_path, 'extents', 'mcaq_surf.tif')
    mc_bot = os.path.join(test_data_path, 'extents', 'mccu_surf.tif')
    lc_top = os.path.join(test_data_path, 'extents', 'lcaq_surf.tif')
    lc_bot = os.path.join(test_data_path, 'extents', 'lccu_surf.tif')
    meras_shp = os.path.join(test_data_path, 'extents', 'MERAS_Extent.shp')

    # make a swuds object and process it
    wu = swuds.Swuds(xlsx=swuds_input, sheet=worksheet, csvfile=outcsv)
    wu.sort_sites()
    wu.reproject()
    wu.apply_footprint(meras_shp, outshp=outputshp)
    wu.assign_missing_elev(dem)

    mc = ['middle_claiborne', mc_top, mc_bot]
    lc = ['lower_claiborne', lc_top, lc_bot]
    wu.make_production_zones([mc, lc])
    wu.assign_monthly_production(processed_csv)

    assert os.path.exists(outputfile)
    assert os.path.exists(outputshp)
    assert os.path.exists(processed_csv)
    
