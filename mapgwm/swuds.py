import os
import sys
import pandas as pd
import numpy as np
from osgeo import ogr
from shapely.geometry import Point
from gisutils import df2shp, project, raster, shp2df
from mapgwm.lookups import aq_codes_dict

class swuds:
    """ Code for preprocessing non-agricultural water use information
    into clean CSV input to MODFLOW setup. Includes logic to fill missing  data,
    as data at many sites are limited to survey years (e.g. 2010 and 2015).
    """

    def __init__(self, xlsx=None, sheet=None, csvfile=None, cols='default'):
        """ Constructor for the swuds class. Class methods will help pre-process
        water-use data for MAP models.  Dataframe produced by constructor
        has original SWUDS data, useful for debugging.

        Parameters
        ----------
        xlsx: str
            Path to xlsx data to be read.  If xlsx file is passed, then a
            selected worksheet from it will be converted to a csv file unless
            the csvfile parameter is specified as None.
        csvfile: str
            Path to csv file with data (if xlsx if None) or to a csvfile
            that is created from the selected worksheet.  If xlsx is None,
            then csvfile must be provided.
        sheet: str
            Name of worksheet in xlsx to be read, ignored if xlsx is None.
        cols: list of str
            List of columns to read from xlsx or csv, if None all columns are read.
            If 'default' (which is the default if nothing is specified) the default
            list of columns coded in the script is read.  

        Attributes
        ----------
        df_swuds: pandas dataframe
            pandas dataframe read from spreadsheet or csv.  Manipulated by other
            methods of the class
        aquifer_names: dict
            dictionary of aquifer names keyed by NWIS codes, read using import 
            statement from mapgwm.lookups
        regional_aquifers: dict
            dictionary of regional aquifers keyed by NWIS codes, read using import
            statement from mapgwm.lookups
        """
        defaultcols=["SITE_NO","WATER_CD","FROM_DEC_LAT_VA","FROM_DEC_LONG_VA","FROM_WELL_DEPTH_VA",
          "FROM_NAT_WATER_USE_CD","FROM_NAT_AQFR_CD","FROM_NAT_AQFR_NM","FROM_AQFR_CD",
          "FROM_AQFR_NM","YEAR","SALINITY_CD","JAN_VAL","FEB_VAL","MAR_VAL","APR_VAL","MAY_VAL",
          "JUN_VAL","JUL_VAL","AUG_VAL","SEP_VAL","OCT_VAL","NOV_VAL","DEC_VAL","ANNUAL_VAL",
          "FROM_STATE_NM","FROM_COUNTY_NM","FROM_CONSTRUCTION_DT","FROM_INVENTORY_DT"]

        if isinstance(cols, list):
            usecols = list
        elif cols is not None:
            usecols = defaultcols
        else:
            usecols = None

        if xlsx is None and csvfile is None:
            sys.exit('both xlsx and csvfile cannont be None for SWUDS object')
        elif xlsx is None:
            self.df_swuds = pd.read_csv(csvfile, usecols=usecols)
        else:
            self.df_swuds = pd.read_excel(xlsx, sheet_name=sheet, usecols=usecols)
            if csvfile is not None:
                self.df_swuds.to_csv(csvfile, index=False)

        self.df_swuds.columns = self.df_swuds.columns.str.upper()
        # remove trailing spaces from code names
        if 'FROM_AQFR_CD' in self.df_swuds.columns:
            self.df_swuds["FROM_AQFR_CD"] = self.df_swuds["FROM_AQFR_CD"].str.strip()

        # make depths floats
        if 'SCREEN_BOT' in self.df_swuds.columns:
            self.df_swuds['SCREEN_BOT'] = pd.to_numeric(self.df_swuds['FROM_WELL_DEPTH_VA'], errors='coerce')
        if 'FROM_ALT_VA' in self.df_swuds.columns:
            self.df_swuds['FROM_ALT_VA'] = pd.to_numeric(self.df_swuds['FROM_ALT_VA'], errors='coerce')

        # other attributes for object
        self.aquifer_names = aq_codes_dict['aquifer_code_names']
        self.regional_aquifers = aq_codes_dict['regional_aquifer']


    def sort_sites(self, secondarysort=None):
        """ Sort the dataframe by site number and quantity, or passed parameter

        Parameters
        ----------
        secondarysort: str
            variable in dataframe to sort SITE_NO groups, default is None
        """
        if secondarysort is not None:
            self.df_swuds = self.df_swuds.sort_values(['SITE_NO', secondarysort]).groupby(['SITE_NO']).last().reset_index()
        else:
            self.df_swuds = self.df_swuds.sort_values(['SITE_NO']).groupby(['SITE_NO']).last().reset_index()

    def reproject(self, to_code=5070):
        """ Reproject from lat/lon to Albers (or passed epsg code) using gisutils

        Parameters
        ----------
        to_code: int
            valid epsg code for coordinate system to be projected into.  Defaults to 5070 - Albers
        """
        x_reprj, y_reprj = project(zip(self.df_swuds['FROM_DEC_LONG_VA'], self.df_swuds['FROM_DEC_LAT_VA']), 
                                 'epsg:4269', 
                                 'epsg:{0}'.format(to_code))
        self.df_swuds['x_{0}'.format(to_code)] = x_reprj
        self.df_swuds['y_{0}'.format(to_code)] = y_reprj
        self.df_swuds['geometry'] = [Point(x, y) for x, y in zip(x_reprj, y_reprj)]

        # drop if no location information
        self.df_swuds.dropna(subset=['x_{0}'.format(to_code), 'y_{0}'.format(to_code)], axis=0, inplace=True)
    

    def apply_footprint(self, bounding_shp, epsg=5070, outshp=None):
        """ Keep sites in the df_swuds pandas dataframe that fall
        into the passed bounding shapefile polygon. Requires
        that df_swuds dataframe has a geometry column as assigned
        in the reproject method. Method uses the first polygon
        in the bounding shapefile if it contains more than one.
        todo: dissolve multiple polygons in bounding_shp?

        Parameters
        ----------
        bounding_shp: str
            path to shapefile with footprint for current analysis
        epsg: int
            valid epsg code for coordinate system for df_swuds.  Defaults to 5070 - Albers
        outshp: str
            optional path to output shapefile with points within the extent
        """
        g2 = shp2df(bounding_shp, dest_crs=epsg)
        poly = g2.geometry.values[0]
        within = [g.within(poly) for g in self.df_swuds['geometry']]
        self.df_swuds = self.df_swuds.loc[within].copy()
        if outshp is not None:
            df2shp(self.df_swuds, outshp)


    # # cull sites to those within the MERAS footprint
    # # meras_bound = '../../../meras3/source_data/extents/meras_extent.shp'
    # delta_bound = '../source_data/extents/delta_active_area_boundary.shp'
    # df_swuds = gpd.GeoDataFrame(df_swuds, crs='epsg:5070')
    # g2 = gpd.GeoDataFrame.from_file(delta_bound)
    # poly = g2.geometry.values[0]
    # within = [g.within(poly) for g in df_swuds.geometry]
    # df_swuds = df_swuds.loc[within].copy()
    # print('national aquifer codes within remaining GW sites: {}'.format(set(df_swuds['FROM_NAT_AQFR_CD'].values.tolist())))
    # # print('number of GW sites in MRVA, Embayment, Alluvial aquivers and NaN: {}'.format(len(df_swuds)))

    # # sample elevations for well with no elevation info
    # no_elev = df_swuds['FROM_ALT_VA'].isnull()
    # x_no_elev = df_swuds.loc[no_elev, 'x_5070'].values
    # y_no_elev = df_swuds.loc[no_elev, 'y_5070'].values
    # elevs = raster.get_values_at_points('../source_data/layers/model_top.tif',
    #                                     x=x_no_elev,
    #                                     y=y_no_elev
    #                                     )
    # df_swuds.loc[no_elev, 'FROM_ALT_VA'] = elevs
    # assert not df_swuds['FROM_ALT_VA'].isnull().any()

    # # make dictionaries of location and elevation information by well
    # locations = dict(list(zip(df_swuds['SITE_NO'], list(zip(df_swuds['x_5070'], df_swuds['y_5070'])))))
    # depths_m = dict(list(zip(df_swuds['SITE_NO'], df_swuds['SCREEN_BOT'] * 0.3048)))
    # well_elevations_m = dict(zip(df_swuds['SITE_NO'], df_swuds['FROM_ALT_VA'] * 0.3048))

    # # get tops and bottoms of estimated production intervals at each well
    # # make dictionaries to lookup by well
    # lcaq_top = '../source_data/water_use/pz_lcaqp_top/pz_lcaqp_top.tif'
    # lcaq_bot = '../source_data/water_use/pz_lcaqp_bot/pz_lcaqp_bot.tif'
    # mcaq_top = '../source_data/water_use/pz_mcaqp_top/pz_mcaqp_top.tif'
    # mcaq_bot = '../source_data/water_use/pz_mcaqp_bot/pz_mcaqp_bot.tif'
    # x_5070 = df_swuds['x_5070'].values
    # y_5070 = df_swuds['y_5070'].values
    # lc_top = raster.get_values_at_points(lcaq_top,
    #                                     x=x_5070,
    #                                     y=y_5070) * 0.3048
    # lc_top_m = dict(zip(df_swuds['SITE_NO'], lc_top))
    # lc_bot = raster.get_values_at_points(lcaq_bot,
    #                                     x=x_5070,
    #                                     y=y_5070) * 0.3048
    # lc_bot_m = dict(zip(df_swuds['SITE_NO'], lc_bot))
    # mc_top = raster.get_values_at_points(mcaq_top,
    #                                     x=x_5070,
    #                                     y=y_5070) * 0.3048
    # mc_top_m = dict(zip(df_swuds['SITE_NO'], mc_top))
    # mc_bot = raster.get_values_at_points(mcaq_bot,
    #                                     x=x_5070,
    #                                     y=y_5070) * 0.3048
    # mc_bot_m = dict(zip(df_swuds['SITE_NO'], mc_bot))

    # monthly_cols = ['JAN_VAL', 'FEB_VAL', 'MAR_VAL', 'APR_VAL', 'MAY_VAL',
    #                 'JUN_VAL', 'JUL_VAL', 'AUG_VAL', 'SEP_VAL', 'OCT_VAL',
    #                 'NOV_VAL', 'DEC_VAL']
    # sim_start_dt = '2008-01-01'
    # sim_end_dt = '2017-12-31'
    # default_screen_len = 50 * 0.3048  # meters
    # print('number of sites from SWUDS: {}'.format(len(df_swuds)))
    # df_swuds = df_swuds.loc[(df_swuds['WATER_CD'] == 'GW') & ~(df_swuds['FROM_NAT_WATER_USE_CD'] == 'IR') &
    #                         ~(df_swuds['FROM_NAT_WATER_USE_CD'] == 'AQ') & ~(df_swuds['FROM_NAT_WATER_USE_CD'] == 'TE')]
    # print('number of GW sites, excluding IR and AQ = {}'.format(len(df_swuds)))
    # df_swuds = df_swuds.loc[~np.isnan(df_swuds['FROM_DEC_LAT_VA'])]
    # print('number of GW sites, excluding IR and AQ, with lat, lon information = {}'.format(len(df_swuds)))
    # df_swuds.head()

    # for c in monthly_cols:
    #     idx = df_swuds.loc[df_swuds[c].isnull()].index.values
    #     df_swuds.loc[idx, c] = df_swuds.loc[idx, 'ANNUAL_VAL']
    # df_swuds.head()

    # # reshape dataframe to have monthly values in same column
    # stacked = pd.DataFrame(df_swuds[monthly_cols].stack())
    # stacked.reset_index(inplace=True)
    # stacked.rename(columns={'level_1': 'month',
    #                         0: 'q_monthly'}, inplace=True)
    # stacked.q_monthly = stacked.q_monthly
    # stacked.index = stacked.level_0
    # stacked = stacked.join(df_swuds)
    # keep_cols = [c for c in stacked.columns if c not in monthly_cols]
    # stacked = stacked[keep_cols]
    # month = {name: i + 1 for i, name in enumerate(monthly_cols)}
    # dates = ['{}-{:02d}'.format(year, month[month_column_name])
    #         for year, month_column_name in zip(stacked.YEAR, stacked.month)]
    # stacked['datetime'] = pd.to_datetime(dates)
    # stacked.sort_values(by=['SITE_NO', 'datetime'], inplace=True)

    # groups = stacked.groupby('SITE_NO')
    # all_groups = []
    # for site_no, group in groups:
    #     group = group.copy()
    #     group.index = pd.to_datetime(group['datetime'])
    #     start_date = pd.Timestamp(sim_start_dt)
    #     end_date = pd.Timestamp(sim_end_dt)

    #     monthly_values_2010 = group.loc[group.datetime.dt.year == 2010]
    #     monthly_values_2010 = dict(zip(monthly_values_2010.datetime.dt.month,
    #                                 monthly_values_2010.q_monthly))
    #     avg_monthly_values = group.groupby(group.index.month).mean().q_monthly.to_dict()
    #     q_mean = group.q_monthly.mean()

    #     # reindex the site data to include all months for simulation period
    #     all_dates = pd.date_range(start_date, end_date, freq='MS')
    #     group = group.reindex(all_dates)
    #     # fill empty dates
    #     q = []
    #     for month, q_monthly in zip(group.index.month, group.q_monthly):
    #         # try to use 2010 values if they exist
    #         if np.isnan(q_monthly):
    #             q_monthly = monthly_values_2010.get(month, np.nan)
    #         # otherwise take the average value for each month
    #         if np.isnan(q_monthly):
    #             q_monthly = avg_monthly_values[month]
    #         # fill missing months with the mean value for the site
    #         if np.isnan(q_monthly):
    #             q_monthly = q_mean
    #         q.append(q_monthly)
    #     # group['q'] = [monthly_values_2010[month] if np.isnan(q_monthly) else q_monthly
    #     #            for month, q_monthly in zip(group.index.month, group.q_monthly)]
    #     group['q'] = q
    #     group['q'] = group['q'] * 3785.4  # convert from mgd to cubic m per d

    #     group['site_no'] = site_no
    #     group['well_elev_m'] = well_elevations_m[site_no]
    #     group['depth_m'] = depths_m[site_no]
    #     well_botm_depth = well_elevations_m[site_no] - depths_m[site_no]
    #     group['x_5070'] = np.nanmin(group['x_5070'])
    #     group['y_5070'] = np.nanmin(group['y_5070'])

    #     lc_top_at_well = lc_top_m[site_no]
    #     lc_bot_at_well = lc_bot_m[site_no]
    #     mc_top_at_well = mc_top_m[site_no]
    #     mc_bot_at_well = mc_bot_m[site_no]

    #     lccu_midpt = np.mean((lc_top_at_well, mc_bot_at_well))  # Lower Claiborne confining unit mid-point
    #     production_zone = 'middle claiborne'
    #     if not np.isnan(lccu_midpt):
    #         if well_botm_depth < lccu_midpt:
    #             production_zone = 'lower claiborne'
    #     elif np.isnan(mc_bot_at_well) and np.isnan(lc_top_at_well):
    #         raise NotImplementedError("Condition of no pumping surface not handled")
    #     elif np.isnan(mc_bot_at_well):
    #         production_zone = 'lower claiborne'
    #     elif np.isnan(lc_bot_at_well):
    #         production_zone = 'middle claiborne'
    #     group['production_zone'] = production_zone

    #     # assign open intervals to wells coded as middle or lower claiborne
    #     # based on estimated production zones
    #     if production_zone == 'middle claiborne':
    #         group['screen_bot_m'] = mc_bot_at_well
    #         group['screen_top_m'] = mc_top_at_well
    #     elif production_zone == 'lower claiborne':
    #         group['screen_bot_m'] = lc_bot_at_well
    #         group['screen_top_m'] = lc_top_at_well

    #     # assign open intervals to wells not coded as middle or lower claiborne
    #     # based on well depth
    #     group['open_int_method'] = 'production zone'
    #     aquifer_name = aquifer_names.get(group["FROM_AQFR_CD"].values[0], 'unnamed')
    #     if aquifer_name not in {'middle claiborne', 'lower claiborne'}:
    #         group['screen_bot_m'] = well_elevations_m[site_no] - depths_m[site_no]
    #         group['screen_top_m'] = well_elevations_m[site_no] - depths_m[site_no] + default_screen_len
    #         group['open_int_method'] = 'well depth'

    #     # add aquifer name
    #     group['aquifer_name'] = aquifer_name

    #     # lc_ovlp = overlaps([lc_bot,lc_top],[group['screen_bot_m'][0],group['screen_top_m'][0]])
    #     # mc_ovlp = overlaps([mc_bot,mc_top],[group['screen_bot_m'][0],group['screen_top_m'][0]])

    #     cols = ['site_no', 'q', 'q_monthly', 'month', 'well_elev_m', 'depth_m',
    #             'screen_bot_m', 'screen_top_m', 'x_5070', 'y_5070']
    #     all_groups.append(group[cols])

    # df = pd.concat(all_groups)
    # df['datetime'] = df.index
    # outfile = '../source_data/water_use/swuds_nonTE_{0}_{1}.csv'.format(sim_start_dt, sim_end_dt)
    # df.to_csv(outfile, index=False)
    # j = 2


    def overlaps(self, a, b):
        """ Return the amount of overlap, in bp
        between a and b.

        Parameters
        ----------
        a: 
        b: 

        Returns
        -------
        int:
            If >0, the number of bp of overlap
            If 0,  they are book-ended.
            If <0, the distance in bp between them
        """

        return min(a[1], b[1]) - max(a[0], b[0])


if __name__ == '__main__':

    home = os.getcwd()
    data_path = os.path.join(os.path.dirname(home), 'working_data')
    swuds_input = os.path.join(data_path, 'LMG-withdrawals-2000-2018.xlsx')
    worksheet = 'LMG-withdrawals-2000-2018'
    outcsv = os.path.join(data_path, 'swuds_nonTE.csv')

    meras_shp = os.path.join(data_path, 'MERAS_Extent.shp')
    wu_shp = os.path.join(data_path, 'WU_points.shp')
    # make a swuds object
    # swuds = swuds(xlsx=swuds_input, sheet=worksheet, csvfile=outcsv)
    swuds = swuds(xlsx=None, sheet=None, csvfile=outcsv, cols=None)
    swuds.sort_sites()
    swuds.reproject()
    swuds.apply_footprint(meras_shp, outshp=wu_shp)
    print(swuds.df_swuds.head())

    # df_swuds = pd.read_csv(outcsv)
    # monthly_cols = ['JAN_VAL', 'FEB_VAL', 'MAR_VAL', 'APR_VAL', 'MAY_VAL',
    #                 'JUN_VAL', 'JUL_VAL', 'AUG_VAL', 'SEP_VAL', 'OCT_VAL',
    #                 'NOV_VAL', 'DEC_VAL']
    # sim_start_dt = '2008-01-01'
    # sim_end_dt = '2017-12-31'