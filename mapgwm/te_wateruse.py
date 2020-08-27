"""
Code for preprocessing thermo-electric water use information
into clean CSV input to MODFLOW setup.

todo: refactor swuds code to functions with numpy docstrings
The docstrings should describe the input and output parameters
and summarize what the function does.

todo: add estimated production intervals from Lynn Torak's surfaces
(A Leaf has started on this in a Notebook)

todo: add test with snippet of source data
See the MODFLOW6 models Wiki in Teams for the source data

"""
import os, sys
import flopy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import geopandas as gpd
from osgeo import ogr
from shapely.geometry import Point
from gisutils import df2shp, project, raster
import yaml

def xlsx2csv(xlsx,sheet,cols,outcsv,header=0,skiprows=None):
    df_swuds = pd.read_excel(xlsx,sheet_name=sheet,header=header,usecols=cols)
    df_swuds.to_csv(outcsv, index=False)
    return(df_swuds)

# convert 2010 thermoelectric .xlsx to .csv (for quicker reading)
xlsx_te = '../source_data/water_use/2010_thermo_model_estimates.xlsx'
outcsv_2010 = '../source_data/water_use/2010_thermo_model_estimates.csv'
cols=['PLANT CODE', 'NAME OF WATER SOURCE',
      'Latitude, decimal degrees','Longitude, decimal degrees',
      'USGS WATER SOURCE CODE','USGS-Estimated Annual Withdrawal (Mgal/d)',
      'USGS-Estimated Annual Minimum Withdrawal (Mgal/d)',
      'USGS-Estimated Annual Maximum Withdrawal (Mgal/d)']
sheet = 'Report_table_UPDATED'
skiprows = [0,1]
header = 2

# convert 2015 thermoelectric .xlsx to .csv (for quicker reading)
xlsx_te = '../source_data/water_use/2015_te_model_estimates_lat.long_comids.xlsx'
outcsv_2015 = '../source_data/water_use/2015_thermo_model_estimates.csv'
cols=['EIA_PLANT_ID', 'NAME_OF_WATER_SOURCE','LATITUDE','LONGITUDE',
      'WATER_SOURCE_CODE','WITHDRAWAL','MIN_WITHDRAWAL','MAX_WITHDRAWAL']
sheet = '2015_ANNUAL_WD_CU'
header = 0

# df_te_2010 = xlsx2csv(xlsx_te,sheet,cols,outcsv_2010,header,skiprows)
# df_te_2015 = xlsx2csv(xlsx_te,sheet,cols,outcsv_2015,header)

df_2010 = pd.read_csv(outcsv_2010)
df_2010.rename(columns={'NAME OF WATER SOURCE':'src_name', 'Latitude, decimal degrees':'lat_dd',
                      'Longitude, decimal degrees':'lon_dd', 'USGS WATER SOURCE CODE':'src_code',
                      'USGS-Estimated Annual Withdrawal (Mgal/d)':'2010_q_mgd',
                      'USGS-Estimated Annual Minimum Withdrawal (Mgal/d)':'2010_min_q_mgd',
                      'USGS-Estimated Annual Maximum Withdrawal (Mgal/d)':'2010_max_q_mgd',
                      'PLANT CODE':'plant_code'},
               inplace=True)

# reproject from lat/lon to albers
x_5070, y_5070 = project(zip(df_2010['lon_dd'], df_2010['lat_dd']), 'epsg:4269', 'epsg:5070')
df_2010['x_5070'] = x_5070
df_2010['y_5070'] = y_5070
df_2010['geometry'] = [Point(x, y) for x, y in zip(x_5070, y_5070)]

# drop wells with no location information (for now)
df_2010.dropna(subset=['x_5070', 'y_5070'], axis=0, inplace=True)

# cull sites to those within the Delta footprint
# meras_bound = '../../../meras3/source_data/extents/meras_extent.shp'
delta_bound = '../source_data/extents/delta_active_area_boundary.shp'
df_2010 = gpd.GeoDataFrame(df_2010, crs='epsg:5070')
g2 = gpd.GeoDataFrame.from_file(delta_bound)
poly = g2.geometry.values[0]
within = [g.within(poly) for g in df_2010.geometry]
df_2010 = df_2010.loc[within].copy()

# cull sites to groundwater sites
df_2010['src_name'] = df_2010['src_name'].str.lower()
df_2010 = df_2010[(df_2010['src_code']=='GW') | (df_2010['src_name'].str.contains('well'))]

df_2015 = pd.read_csv(outcsv_2015)
df_2015.rename(columns={'NAME_OF_WATER_SOURCE':'src_name', 'LATITUDE':'lat_dd',
                        'LONGITUDE':'lon_dd', 'WATER_SOURCE_CODE':'src_code',
                        'WITHDRAWAL':'2015_q_mgd', 'MIN_WITHDRAWAL':'2015_min_q_mgd',
                        'MAX_WITHDRAWAL':'2015_max_q_mgd', 'EIA_PLANT_ID':'plant_code'},
               inplace=True)

# reproject from lat/lon to albers
x_5070, y_5070 = project(zip(df_2015['lon_dd'], df_2015['lat_dd']), 'epsg:4269', 'epsg:5070')
df_2015['x_5070'] = x_5070
df_2015['y_5070'] = y_5070
df_2015['geometry'] = [Point(x, y) for x, y in zip(x_5070, y_5070)]

# drop wells with no location information (for now)
df_2015.dropna(subset=['x_5070', 'y_5070'], axis=0, inplace=True)

# cull sites to those within the Delta footprint
df_2015 = gpd.GeoDataFrame(df_2015, crs='epsg:5070')
poly = g2.geometry.values[0]
within = [g.within(poly) for g in df_2015.geometry]
df_2015 = df_2015.loc[within].copy()

# cull sites to groundwater sites
df_2015['src_name'] = df_2015['src_name'].str.lower()
df_2015 = df_2015[(df_2015['src_code']=='GW') | (df_2015['src_name'].str.contains('well'))]

# merge 2015 annual withdrawals with 2010 dataframe
df = pd.merge(df_2010, df_2015[['plant_code','2015_q_mgd']], left_on=['plant_code'],
              right_on=['plant_code'], how='left')

# get tops and bottoms of estimated production intervals at each well
# make dictionaries to lookup by well
mcaq_top = '../source_data/water_use/pz_mcaqp_top/pz_mcaqp_top.tif'
mcaq_bot = '../source_data/water_use/pz_mcaqp_bot/pz_mcaqp_bot.tif'
x_5070 = df['x_5070'].values
y_5070 = df['y_5070'].values
mc_top = raster.get_values_at_points(mcaq_top,
                                     x=x_5070,
                                     y=y_5070)*0.3048
mc_top_m = dict(zip(df['plant_code'], mc_top))
mc_bot = raster.get_values_at_points(mcaq_bot,
                                     x=x_5070,
                                     y=y_5070)*0.3048
mc_bot_m = dict(zip(df['plant_code'], mc_bot))

sim_start_dt = '2008-01-01'
sim_end_dt = '2017-12-31'
df['start_date'] = pd.to_datetime('1/1/2010')
groups = df.groupby('plant_code')
all_groups = []
for plant_code, group in groups:
    group = group.copy()
    group.index = group['start_date']
    start_date = pd.Timestamp(sim_start_dt)
    end_date = pd.Timestamp(sim_end_dt)

    monthly_values_2010 = group.loc[group.start_date == pd.to_datetime('1/1/2010'),'2010_q_mgd'].values[0]
    monthly_values_2015 = group.loc[group.start_date == pd.to_datetime('1/1/2010'),'2015_q_mgd'].values[0]

    all_dates = pd.date_range(start_date, end_date, freq='MS')
    group = group.reindex(all_dates)

    # fill empty dates
    mask = (group.index < pd.to_datetime('1/1/2013'))
    group.loc[mask,'q_m3d'] = monthly_values_2010*3785.4
    mask = (group.index > pd.to_datetime('12/31/2012'))
    group.loc[mask,'q_m3d'] = monthly_values_2015*3785.4

    mc_top_at_well = mc_top_m[plant_code]
    mc_bot_at_well = mc_bot_m[plant_code]
    group['screen_top_m'] = mc_top_at_well
    group['screen_bot_m'] = mc_bot_at_well
    group['x_5070'] = np.nanmin(group['x_5070'])
    group['y_5070'] = np.nanmin(group['y_5070'])
    group['plant_code'] = plant_code

    cols = ['plant_code', 'q_m3d', 'screen_bot_m',
            'screen_top_m', 'x_5070', 'y_5070']
    all_groups.append(group[cols])

df_all = pd.concat(all_groups)
df_all['datetime'] = df_all.index
outfile = '../source_data/water_use/TE_model_est_{0}_{1}.csv'.format(sim_start_dt,sim_end_dt)
df_all.to_csv(outfile, index=False)
j=2