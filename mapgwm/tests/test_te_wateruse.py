import os
from mapgwm import te_wateruse


def test_preprocess_thermal(test_output_folder, test_data_path):
   
    # output
    processed_csv = os.path.join(test_output_folder, 'processed_thermal.csv')
    wu_shp = os.path.join(test_output_folder, 'TE_points.shp')
    outcsv_2010 = os.path.join(test_output_folder, '2010_thermo_model_estimates.csv')
    outcsv_2015 = os.path.join(test_output_folder, '2015_thermo_model_estimates.csv')
    combined = os.path.join(test_output_folder, 'TE_combined.csv')

    # required data
    xlsx_2010 = os.path.join(test_data_path, 'swuds', '2010_thermo_model_estimates.xlsx')
    xlsx_2015 = os.path.join(test_data_path, 'swuds','2015_te_model_estimates_lat.long_comids.xlsx')
    skiprows = [2, 0]
    worksheets = ['Report_table_UPDATED', '2015_ANNUAL_WD_CU']

    # input rasters and shapefile
    mc_top = os.path.join(test_data_path, 'rasters', 'mcaq_surf.tif')
    mc_bot = os.path.join(test_data_path, 'rasters', 'mccu_surf.tif')
    meras_shp = os.path.join(test_data_path, 'extents', 'MERAS_Extent.shp')

    # make a thermal object
    te = te_wateruse.ThermalUse(xlsx=[xlsx_2010, xlsx_2015], 
                    sheet=worksheets, 
                    csvfile=[outcsv_2010, outcsv_2015], 
                    cols=['2010', '2015'], 
                    skiprows=skiprows, 
                    newheaders=['2010', '2015'], 
                    epsg=5070)

    # write the combined dataframe 
    te.df_swuds.to_csv(combined)

    # processing steps
    te.sort_sites(primarysort='plant_code')
    te.reproject(long='lon_dd', lat='lat_dd', key='plant_code')

    te.apply_footprint(meras_shp, outshp=wu_shp)
    
    mc = ['middle_claiborne', mc_top, mc_bot]
    te.make_production_zones([mc], key='plant_code')
    
    te.assign_monthly_production(processed_csv, year_cols=['2010_q_mgd', '2015_q_mgd'])

    #check if output are written
    assert os.path.exists(outcsv_2010)
    assert os.path.exists(outcsv_2015)
    assert os.path.exists(combined)
    assert os.path.exists(wu_shp)
    assert os.path.exists(processed_csv)
    
    
