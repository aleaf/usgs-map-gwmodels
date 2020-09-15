from collections import defaultdict
import os
import re
import sys

import numpy as np
import pandas as pd
from shapely.geometry import Point
import yaml

from gisutils import df2shp, project, raster, shp2df
from mapgwm.lookups import aq_codes_dict
from mapgwm.swuds import Swuds


class ThermalUse(Swuds):
    """ Code for preprocessing thermal electric (TE) water use information
    into clean CSV input to MODFLOW setup. Original functions were hard-coded
    for 2010 and 2015 data, and the input data from these years had different
    column headings.  This class repeats that code, but it should be general
    enough to use other data.

    This is a derived class ecause many of the functions needed to pre-process
    the data are in the Swuds() class.
    """

    def __init__(self, xlsx=None, sheet=None, csvfile=None, cols=['2010', '2015'], 
                 skiprows=None, newheaders=['2010', '2015'], keys=['plant_code', 'plant_code'],
                 epsg=5070):
        """ Constructor for the ThermalUse class. Unlike the Swuds() class,
        ThermalUse will accept a list of xlsx names, sheetnames, and csvfiles.
        If multiple input datasets are sent, they are read in and combined to a single
        dataframe. The cols and header parameters map the columns read to common 
        names used for the join.

        Parameters
        ----------
        xlsx: str
            Path to xlsx data to be read.  If xlsx file is passed, then a
            selected worksheet from it will be converted to a csv file unless
            the csvfile parameter is specified as None.  If a list is passed,
            each will be read in, converted, and the dataframes combined.
        csvfile: str
            Path to csv file with data (if xlsx if None) or to a csvfile
            that is created from the selected worksheet.  If xlsx is None,
            then csvfile must be provided. If a list is passed,
            each will be read in or used for output of the xlsx. Individual
            dataframes will be combined.
        sheet: str
            Name of worksheet(s) in xlsx to be read, ignored if xlsx is None.
            Must be a list of xlsx is a list.
        cols: list of str
            List of columns to read from xlsx or csv, if None all columns are read.
            If '['2010', '2015']' (which is the default if nothing is specified) the default
            list of columns for 2010 and 2015 hard-coded in the script are used.  This
            parameter will be a list of lists of strings if xlsx is a list and defaults or
            None is not used.
        skiprows: int
            Number of rows to skip when reading in the xlsx file. Must be a list of ints
            if xlsx is a list.
        newheaders: dict
            Dictionary converting column names in input to common columns 
            used to combine data frames. If the default is used, 
            dictionaries for 2010 and 2015 hard-coded 
            in the script are used.  This parameter will be a list of
            dictionaries if xlsx is a list and the default is not used. 
            If None is passed the column names are not changed.
        keys: str
            keys used to join dataframes into a single dataframe, pass None if only
            one dataset is read in.  Defaults to plant_code used in 2010, 2015 datasets.
        epsg: int
            valid epsg code for coordinate system used in column name in df_swuds.  
            Defaults to 5070 - Albers

        Attributes
        ----------
        same as the Swuds() class.

        """

        # set some attributes
        self.monthly_cols = ['JAN_VAL', 'FEB_VAL', 'MAR_VAL', 'APR_VAL', 'MAY_VAL',
                             'JUN_VAL', 'JUL_VAL', 'AUG_VAL', 'SEP_VAL', 'OCT_VAL',
                             'NOV_VAL', 'DEC_VAL']
        self.aquifer_names = aq_codes_dict['aquifer_code_names']
        self.regional_aquifers = aq_codes_dict['regional_aquifer']
        self.sim_start_dt = '2008-01-01'
        self.sim_end_dt = '2017-12-31'
        self.default_screen_len = 50. * 0.3028
        self.epsg = epsg
        self.prod_zone_top_m = defaultdict(dict)
        self.prod_zone_bot_m = defaultdict(dict)
        self.locations = dict()

        # set defaults for 2010 and 2015
        default2010 = ['PLANT CODE', 'NAME OF WATER SOURCE',
                       'Latitude, decimal degrees','Longitude, decimal degrees',
                       'USGS WATER SOURCE CODE','USGS-Estimated Annual Withdrawal (Mgal/d)',
                       'USGS-Estimated Annual Minimum Withdrawal (Mgal/d)',
                       'USGS-Estimated Annual Maximum Withdrawal (Mgal/d)']
        rename2010 = {'NAME OF WATER SOURCE':'src_name', 'Latitude, decimal degrees':'lat_dd',
                      'Longitude, decimal degrees':'lon_dd', 'USGS WATER SOURCE CODE':'src_code',
                      'USGS-Estimated Annual Withdrawal (Mgal/d)':'2010_q_mgd',
                      'USGS-Estimated Annual Minimum Withdrawal (Mgal/d)':'2010_min_q_mgd',
                      'USGS-Estimated Annual Maximum Withdrawal (Mgal/d)':'2010_max_q_mgd',
                      'PLANT CODE':'plant_code'}
        default2015 = ['EIA_PLANT_ID', 'NAME_OF_WATER_SOURCE','LATITUDE','LONGITUDE',
                       'WATER_SOURCE_CODE','WITHDRAWAL','MIN_WITHDRAWAL','MAX_WITHDRAWAL']
        rename2015 = {'NAME_OF_WATER_SOURCE':'src_name', 'LATITUDE':'lat_dd',
                      'LONGITUDE':'lon_dd', 'WATER_SOURCE_CODE':'src_code',
                      'WITHDRAWAL':'2015_q_mgd', 'MIN_WITHDRAWAL':'2015_min_q_mgd',
                      'MAX_WITHDRAWAL':'2015_max_q_mgd', 'EIA_PLANT_ID':'plant_code'}

        if not isinstance(xlsx, list):
            xlsx = [xlsx]
            csvfile = [csvfile]
            sheet = [sheet]
            cols = [cols]
            newheaders = [newheaders]
            skiprows = [skiprows]
            keys = [keys]

        tes = []

        for i in range(0, len(xlsx)):
            if isinstance(cols[i], list):
                usecols = cols[i]
            elif re.match('2010', cols[i]):
                usecols = default2010
            elif re.match('2015', cols[i]):
                usecols = default2015
            else:
                usecols = None

            if xlsx[i] is None and csvfile[i] is None:
                sys.exit('both xlsx and csvfile cannot be None for ThermalUse object')
            elif xlsx is None:
                te = pd.read_csv(csvfile[i], usecols=usecols)
            else:
                te = pd.read_excel(xlsx[i], sheet_name=sheet[i], skiprows=skiprows[i], usecols=usecols)
                if csvfile is not None:
                    te.to_csv(csvfile[i], index=False)

            # te.columns = te.columns.str.upper()
            if isinstance(newheaders[i], dict):
                renamedict = newheaders[i]
            elif re.match('2010', newheaders[i]):
                renamedict = rename2010
            elif re.match('2015', newheaders[i]):
                renamedict = rename2015
            else:
                renamedict = None
            
            if renamedict is not None:
                te.rename(columns=renamedict, inplace=True)

            if keys[i] is not None:
                te.set_index(keys[i], inplace=True)
            tes.append(te)
        self.df_swuds = pd.concat(tes, axis=1)
        self.df_swuds = self.df_swuds.loc[:, ~self.df_swuds.columns.duplicated()]  # get rid of repeated columns from concat



# # convert 2010 thermoelectric .xlsx to .csv (for quicker reading)
# xlsx_te = '../source_data/water_use/2010_thermo_model_estimates.xlsx'
# outcsv_2010 = '../source_data/water_use/2010_thermo_model_estimates.csv'
# cols=['PLANT CODE', 'NAME OF WATER SOURCE',
#       'Latitude, decimal degrees','Longitude, decimal degrees',
#       'USGS WATER SOURCE CODE','USGS-Estimated Annual Withdrawal (Mgal/d)',
#       'USGS-Estimated Annual Minimum Withdrawal (Mgal/d)',
#       'USGS-Estimated Annual Maximum Withdrawal (Mgal/d)']
# sheet = 'Report_table_UPDATED'
# skiprows = [0,1]
# header = 2

# # convert 2015 thermoelectric .xlsx to .csv (for quicker reading)
# xlsx_te = '../source_data/water_use/2015_te_model_estimates_lat.long_comids.xlsx'
# outcsv_2015 = '../source_data/water_use/2015_thermo_model_estimates.csv'
# cols=['EIA_PLANT_ID', 'NAME_OF_WATER_SOURCE','LATITUDE','LONGITUDE',
#       'WATER_SOURCE_CODE','WITHDRAWAL','MIN_WITHDRAWAL','MAX_WITHDRAWAL']
# sheet = '2015_ANNUAL_WD_CU'
# header = 0

# # df_te_2010 = xlsx2csv(xlsx_te,sheet,cols,outcsv_2010,header,skiprows)
# # df_te_2015 = xlsx2csv(xlsx_te,sheet,cols,outcsv_2015,header)

# df_2010 = pd.read_csv(outcsv_2010)
# df_2010.rename(columns={'NAME OF WATER SOURCE':'src_name', 'Latitude, decimal degrees':'lat_dd',
#                       'Longitude, decimal degrees':'lon_dd', 'USGS WATER SOURCE CODE':'src_code',
#                       'USGS-Estimated Annual Withdrawal (Mgal/d)':'2010_q_mgd',
#                       'USGS-Estimated Annual Minimum Withdrawal (Mgal/d)':'2010_min_q_mgd',
#                       'USGS-Estimated Annual Maximum Withdrawal (Mgal/d)':'2010_max_q_mgd',
#                       'PLANT CODE':'plant_code'},
#                inplace=True)

# # reproject from lat/lon to albers
# x_5070, y_5070 = project(zip(df_2010['lon_dd'], df_2010['lat_dd']), 'epsg:4269', 'epsg:5070')
# df_2010['x_5070'] = x_5070
# df_2010['y_5070'] = y_5070
# df_2010['geometry'] = [Point(x, y) for x, y in zip(x_5070, y_5070)]

# # drop wells with no location information (for now)
# df_2010.dropna(subset=['x_5070', 'y_5070'], axis=0, inplace=True)

# # cull sites to those within the Delta footprint
# # meras_bound = '../../../meras3/source_data/extents/meras_extent.shp'
# delta_bound = '../source_data/extents/delta_active_area_boundary.shp'
# df_2010 = gpd.GeoDataFrame(df_2010, crs='epsg:5070')
# g2 = gpd.GeoDataFrame.from_file(delta_bound)
# poly = g2.geometry.values[0]
# within = [g.within(poly) for g in df_2010.geometry]
# df_2010 = df_2010.loc[within].copy()

# # cull sites to groundwater sites
# df_2010['src_name'] = df_2010['src_name'].str.lower()
# df_2010 = df_2010[(df_2010['src_code']=='GW') | (df_2010['src_name'].str.contains('well'))]

# df_2015 = pd.read_csv(outcsv_2015)
# df_2015.rename(columns={'NAME_OF_WATER_SOURCE':'src_name', 'LATITUDE':'lat_dd',
#                         'LONGITUDE':'lon_dd', 'WATER_SOURCE_CODE':'src_code',
#                         'WITHDRAWAL':'2015_q_mgd', 'MIN_WITHDRAWAL':'2015_min_q_mgd',
#                         'MAX_WITHDRAWAL':'2015_max_q_mgd', 'EIA_PLANT_ID':'plant_code'},
#                inplace=True)

# # reproject from lat/lon to albers
# x_5070, y_5070 = project(zip(df_2015['lon_dd'], df_2015['lat_dd']), 'epsg:4269', 'epsg:5070')
# df_2015['x_5070'] = x_5070
# df_2015['y_5070'] = y_5070
# df_2015['geometry'] = [Point(x, y) for x, y in zip(x_5070, y_5070)]

# # drop wells with no location information (for now)
# df_2015.dropna(subset=['x_5070', 'y_5070'], axis=0, inplace=True)

# # cull sites to those within the Delta footprint
# df_2015 = gpd.GeoDataFrame(df_2015, crs='epsg:5070')
# poly = g2.geometry.values[0]
# within = [g.within(poly) for g in df_2015.geometry]
# df_2015 = df_2015.loc[within].copy()

# # cull sites to groundwater sites
# df_2015['src_name'] = df_2015['src_name'].str.lower()
# df_2015 = df_2015[(df_2015['src_code']=='GW') | (df_2015['src_name'].str.contains('well'))]

# # merge 2015 annual withdrawals with 2010 dataframe
# df = pd.merge(df_2010, df_2015[['plant_code','2015_q_mgd']], left_on=['plant_code'],
#               right_on=['plant_code'], how='left')

# # get tops and bottoms of estimated production intervals at each well
# # make dictionaries to lookup by well
# mcaq_top = '../source_data/water_use/pz_mcaqp_top/pz_mcaqp_top.tif'
# mcaq_bot = '../source_data/water_use/pz_mcaqp_bot/pz_mcaqp_bot.tif'
# x_5070 = df['x_5070'].values
# y_5070 = df['y_5070'].values
# mc_top = raster.get_values_at_points(mcaq_top,
#                                      x=x_5070,
#                                      y=y_5070)*0.3048
# mc_top_m = dict(zip(df['plant_code'], mc_top))
# mc_bot = raster.get_values_at_points(mcaq_bot,
#                                      x=x_5070,
#                                      y=y_5070)*0.3048
# mc_bot_m = dict(zip(df['plant_code'], mc_bot))

# sim_start_dt = '2008-01-01'
# sim_end_dt = '2017-12-31'
# df['start_date'] = pd.to_datetime('1/1/2010')
# groups = df.groupby('plant_code')
# all_groups = []
# for plant_code, group in groups:
#     group = group.copy()
#     group.index = group['start_date']
#     start_date = pd.Timestamp(sim_start_dt)
#     end_date = pd.Timestamp(sim_end_dt)

#     monthly_values_2010 = group.loc[group.start_date == pd.to_datetime('1/1/2010'),'2010_q_mgd'].values[0]
#     monthly_values_2015 = group.loc[group.start_date == pd.to_datetime('1/1/2010'),'2015_q_mgd'].values[0]

#     all_dates = pd.date_range(start_date, end_date, freq='MS')
#     group = group.reindex(all_dates)

#     # fill empty dates
#     mask = (group.index < pd.to_datetime('1/1/2013'))
#     group.loc[mask,'q_m3d'] = monthly_values_2010*3785.4
#     mask = (group.index > pd.to_datetime('12/31/2012'))
#     group.loc[mask,'q_m3d'] = monthly_values_2015*3785.4

#     mc_top_at_well = mc_top_m[plant_code]
#     mc_bot_at_well = mc_bot_m[plant_code]
#     group['screen_top_m'] = mc_top_at_well
#     group['screen_bot_m'] = mc_bot_at_well
#     group['x_5070'] = np.nanmin(group['x_5070'])
#     group['y_5070'] = np.nanmin(group['y_5070'])
#     group['plant_code'] = plant_code

#     cols = ['plant_code', 'q_m3d', 'screen_bot_m',
#             'screen_top_m', 'x_5070', 'y_5070']
#     all_groups.append(group[cols])

# df_all = pd.concat(all_groups)
# df_all['datetime'] = df_all.index
# outfile = '../source_data/water_use/TE_model_est_{0}_{1}.csv'.format(sim_start_dt,sim_end_dt)
# df_all.to_csv(outfile, index=False)
# j=2


if __name__ == '__main__':

    home = os.getcwd()
    data_path = os.path.join(os.path.dirname(home), 'working_data')
    xlsx_2010 = os.path.join(data_path, '2010_thermo_model_estimates.xlsx')
    xlsx_2015 = os.path.join(data_path, '2015_te_model_estimates_lat.long_comids.xlsx')
    outcsv_2010 = os.path.join(data_path, '2010_thermo_model_estimates.csv')
    outcsv_2015 = os.path.join(data_path, '2015_thermo_model_estimates.csv')
    skiprows = [2, 0]
    worksheets = ['Report_table_UPDATED', '2015_ANNUAL_WD_CU']

    te = ThermalUse(xlsx=[xlsx_2010, xlsx_2015], 
                    sheet=worksheets, 
                    csvfile=[outcsv_2010, outcsv_2015], 
                    cols=['2010', '2015'], 
                    skiprows=skiprows, 
                    newheaders=['2010', '2015'], 
                    epsg=5070)

    
    dem = os.path.join(data_path, 'dem_mean_elevs.tif')
    mc_top = os.path.join(data_path, 'mcaq_surf.tif')
    mc_bot = os.path.join(data_path, 'mccu_surf.tif')
    lc_top = os.path.join(data_path, 'lcaq_surf.tif')
    lc_bot = os.path.join(data_path, 'lccu_surf.tif')

    meras_shp = os.path.join(data_path, 'MERAS_Extent.shp')
    wu_shp = os.path.join(data_path, 'TE_points.shp')
    
    te.sort_sites(primarysort='plant_code')
    te.reproject(long='lon_dd', lat='lat_dd', key='plant_code')

    te.apply_footprint(meras_shp, outshp=wu_shp)
    print(te.df_swuds.head())
    # swuds.assign_missing_elev(dem)

    # mc = ['middle_claiborne', mc_top, mc_bot]
    # lc = ['lower_claiborne', lc_top, lc_bot]
    # swuds.make_production_zones([mc, lc])
    # swuds.assign_monthly_production(os.path.join(data_path, 'processed_swuds.csv'))
    # print(swuds.df_swuds.head())

    # yml_file = os.path.join(data_path, 'swuds_input.yml')
    # wu = Swuds.from_yaml(yml_file)