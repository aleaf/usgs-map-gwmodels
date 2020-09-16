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


class Swuds:
    """ Code for preprocessing non-agricultural water use information
    into clean CSV input to MODFLOW setup. Class excludes AQ, IR, and TE
    water use from a swuds excel dataset. Includes logic to fill missing  data,
    as data at many sites are limited to survey years (e.g. 2010 and 2015).
    """

    def __init__(self, xlsx=None, sheet=None, csvfile=None, cols='default', epsg=5070):
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
        epsg: int
            valid epsg code for coordinate system used in column name in df_swuds.  
            Defaults to 5070 - Albers

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
        epsg: int
            epsg code used in projections
        monthly_cols: list of str
            list of column names for monthly values
        sim_start_dt: str
            start time for simulation as string 'yyyy-mm-dd'
        sim_end_dt: str
            end time for simulation as string 'yyyy-mm-dd'
        default_screen_len: float
            default screen length in meters
        locations: dict
            dictionary of x,y locations keyed by SITE_NO, added in reproject method
        depths_m: dict
            dictionary of depth keyed by SITE_NO
        well_elevations_m: dict
            dictionary of well elevations keyed by SITE_NO
        prod_zone_top_m: defaultdict(dict)
            defaultdict with production zone top (in meters) for each well
            First key is the production zone name, and second is the SITE_NO
            For example  self.prod_zone_top_m['lower_claiborne']['WEL001'] = top_elev
        prod_zone_bot_m: defaultdict(dict)
            defaultdict with production zone bottom (in meters) for each well
            First key is the production zone name, and second is the SITE_NO
            For example  self.prod_zone_bot_m['lower_claiborne']['WEL001'] = bot_elev
        
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

        # now read in excel or csv file
        defaultcols=["SITE_NO","WATER_CD","FROM_DEC_LAT_VA","FROM_DEC_LONG_VA","FROM_WELL_DEPTH_VA",
          "FROM_ALT_VA","FROM_NAT_WATER_USE_CD","FROM_NAT_AQFR_CD","FROM_NAT_AQFR_NM","FROM_AQFR_CD",
          "FROM_AQFR_NM","YEAR","SALINITY_CD","JAN_VAL","FEB_VAL","MAR_VAL","APR_VAL","MAY_VAL",
          "JUN_VAL","JUL_VAL","AUG_VAL","SEP_VAL","OCT_VAL","NOV_VAL","DEC_VAL","ANNUAL_VAL",
          "FROM_STATE_NM","FROM_COUNTY_NM","FROM_CONSTRUCTION_DT","FROM_INVENTORY_DT"]

        if isinstance(cols, list):
            usecols = cols
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
        if 'FROM_WELL_DEPTH_VA' in self.df_swuds.columns:
            self.df_swuds['SCREEN_BOT'] = pd.to_numeric(self.df_swuds['FROM_WELL_DEPTH_VA'], errors='coerce')
        else:
            self.df_swuds['SCREEN_BOT'] = np.nan
        if 'FROM_ALT_VA' in self.df_swuds.columns:
            self.df_swuds['FROM_ALT_VA'] = pd.to_numeric(self.df_swuds['FROM_ALT_VA'], errors='coerce')

        # make dictionaries
        self.depths_m = dict(list(zip(self.df_swuds['SITE_NO'], self.df_swuds['SCREEN_BOT'] * 0.3048)))
        self.well_elevations_m = dict(zip(self.df_swuds['SITE_NO'], self.df_swuds['FROM_ALT_VA'] * 0.3048))

        
    def sort_sites(self, primarysort='SITE_NO', secondarysort=None):
        """ Sort the dataframe by site number and quantity, or passed parameter

        Parameters
        ----------
        primarysort: str
            string to sort group groups, defaults to SITE_NO
        secondarysort: str
            variable in dataframe to sort SITE_NO groups, default is None
        """
        if secondarysort is not None:
            self.df_swuds = self.df_swuds.sort_values([primarysort, secondarysort]).groupby([primarysort]).last().reset_index()
        else:
            self.df_swuds = self.df_swuds.sort_values([primarysort]).groupby([primarysort]).last().reset_index()

    def reproject(self, long='FROM_DEC_LONG_VA', lat='FROM_DEC_LAT_VA', key='SITE_NO'):
        """ Reproject from lat/lon to self.epsg using gisutils

        Parameters
        ----------
        long: str
            string variable name (column in dataframe) with longitude for input
        lat: str
            string variable name (column in dataframe) with latitude 
        key: str
            key for the dictionary made, defaults to SITE_NO
        """
        x_reprj, y_reprj = project(zip(self.df_swuds[long], self.df_swuds[lat]), 
                                 'epsg:4269', 
                                 'epsg:{0}'.format(self.epsg))
        self.df_swuds['x_{0}'.format(self.epsg)] = x_reprj
        self.df_swuds['y_{0}'.format(self.epsg)] = y_reprj
        self.df_swuds['geometry'] = [Point(x, y) for x, y in zip(x_reprj, y_reprj)]

        # drop entries if no location information
        self.df_swuds.dropna(subset=['x_{0}'.format(self.epsg), 'y_{0}'.format(self.epsg)], axis=0, inplace=True)

        # make dictionary of location 
        self.locations = dict(list(zip(self.df_swuds[key], list(zip(self.df_swuds['x_{0}'.format(self.epsg)], self.df_swuds['y_{0}'.format(self.epsg)])))))
        
    
    def apply_footprint(self, bounding_shp, outshp=None):
        """ Keep sites in the df_swuds pandas dataframe that fall
        into the passed bounding shapefile polygon. Requires
        that df_swuds dataframe has a Point geometry column as assigned
        in the reproject method. Method uses the first polygon
        in the bounding shapefile if it contains more than one.
        todo: dissolve multiple polygons in bounding_shp?

        Parameters
        ----------
        bounding_shp: str
            path to shapefile with footprint for current analysis
        outshp: str
            optional path to output shapefile with points within the footprint
        """
        g2 = shp2df(bounding_shp, dest_crs=self.epsg)
        poly = g2.geometry.values[0]
        within = [g.within(poly) for g in self.df_swuds['geometry']]   # shapely within method
        self.df_swuds = self.df_swuds.loc[within].copy()
        if outshp is not None:
            df2shp(self.df_swuds, outshp)


    def assign_missing_elev(self, top_raster, elev_field='FROM_ALT_VA'):
        """ Use the top of model raster, or land-surface raster,
        to assign the elevation for points where elevation is missing.
        
        Parameters
        ----------
        top_raster: str
            path to raster data set with land surface or model top elevation, used
            to assign missing values to water-use points
        elev_field: str
            field in df_swuds with elevation data, default is 'FROM_ALT_VA'
        """

        no_elev = self.df_swuds['FROM_ALT_VA'].isnull()
        x_no_elev = self.df_swuds.loc[no_elev, 'x_{0}'.format(self.epsg)].values
        y_no_elev = self.df_swuds.loc[no_elev, 'y_{0}'.format(self.epsg)].values
        elevs = raster.get_values_at_points(top_raster,
                                            x=x_no_elev,
                                            y=y_no_elev
                                            )
        self.df_swuds.loc[no_elev, 'FROM_ALT_VA'] = elevs
        assert not self.df_swuds['FROM_ALT_VA'].isnull().any()

    
    def make_production_zones(self, zonelist, key='SITE_NO'):
        """ Make dictionary attributes for production zones.
        These are used to assign individual wells to production zones.
        The defaultdict is keyed by zone_name and then SITE_NO.

        Parameters
        ----------
        zonelist: list of lists
            List of production zone information, each zone requires a 
            list with [zone_name, zone_top, zone_bot]
        zone_name: str
            name assigned to prodcuction zone
        zone_top: str
            path to raster with top of zone
        zone_bot: str
            path to raster to bottom of zone
        key: str
            key (column name) to use in the resulting
            parameter zone dictionaries.  Defaults to SITE_NO
        """

        # if only one list is passed, put it into a list.
        if isinstance(zonelist[0], str):
            zonelist = [zonelist]

        # get tops and bottoms of estimated production intervals at each well
        # make dictionaries to lookup by well
        for z in zonelist:
            name = z[0]
            top_raster = z[1]
            bot_raster = z[2]
            x = self.df_swuds['x_{0}'.format(self.epsg)].values
            y = self.df_swuds['y_{0}'.format(self.epsg)].values
            lc_top = raster.get_values_at_points(top_raster,
                                                x=x,
                                                y=y) * 0.3048
            self.prod_zone_top_m[name] = dict(zip(self.df_swuds[key], lc_top))
            lc_bot = raster.get_values_at_points(bot_raster,
                                                x=x,
                                                y=y) * 0.3048
            self.prod_zone_bot_m[name] = dict(zip(self.df_swuds[key], lc_bot))
            

    def assign_monthly_production(self, outfile):
        """ Assign production wells for water use, skipping IR (irrigation) and
        TE (thermal electric) to production zones.  If production zones are not
        assigned or if the well bottom doesn't fall into a production zone, then
        the screen_top and screen_bot are assigned using well_depth and the
        default screen length.

        Production is given in cubic m per day.
        todo:  add unit conversion parameter so other units can be used?

        Parameters
        ----------
        outfile: str
            path to final processed monthly water-use file with production zone
            information 
        """

        # fill in missing monthly values with annual value
        for c in self.monthly_cols:
            idx = self.df_swuds.loc[self.df_swuds[c].isnull()].index.values
            self.df_swuds.loc[idx, c] = self.df_swuds.loc[idx, 'ANNUAL_VAL']

        # pull out groundwater sites that are not IR, AQ or TE
        self.df_swuds = self.df_swuds.loc[
                                (self.df_swuds['WATER_CD'] == 'GW') & 
                                ~(self.df_swuds['FROM_NAT_WATER_USE_CD'] == 'IR') & 
                                ~(self.df_swuds['FROM_NAT_WATER_USE_CD'] == 'AQ') & 
                                ~(self.df_swuds['FROM_NAT_WATER_USE_CD'] == 'TE')
                                ]


        # reshape dataframe to have monthly values in same column
        stacked = pd.DataFrame(self.df_swuds[self.monthly_cols].stack())
        stacked.reset_index(inplace=True)
        stacked.rename(columns={'level_1': 'month',
                                0: 'q_monthly'}, inplace=True)
        stacked.q_monthly = stacked.q_monthly
        stacked.index = stacked.level_0
        stacked = stacked.join(self.df_swuds)
        keep_cols = [c for c in stacked.columns if c not in self.monthly_cols]
        stacked = stacked[keep_cols]
        month = {name: i + 1 for i, name in enumerate(self.monthly_cols)}
        dates = ['{}-{:02d}'.format(year, month[month_column_name])
                for year, month_column_name in zip(stacked.YEAR, stacked.month)]
        stacked['datetime'] = pd.to_datetime(dates)
        stacked.sort_values(by=['SITE_NO', 'datetime'], inplace=True)

        groups = stacked.groupby('SITE_NO')
        all_groups = []
        for site_no, group in groups:
            group = group.copy()
            group.index = pd.to_datetime(group['datetime'])
            start_date = pd.Timestamp(self.sim_start_dt)
            end_date = pd.Timestamp(self.sim_end_dt)

            monthly_values_2010 = group.loc[group.datetime.dt.year == 2010]
            monthly_values_2010 = dict(zip(monthly_values_2010.datetime.dt.month,
                                        monthly_values_2010.q_monthly))
            avg_monthly_values = group.groupby(group.index.month).mean().q_monthly.to_dict()
            q_mean = group.q_monthly.mean()

            # reindex the site data to include all months for simulation period
            all_dates = pd.date_range(start_date, end_date, freq='MS')
            group = group.reindex(all_dates)
            # fill empty dates
            q = []
            for month, q_monthly in zip(group.index.month, group.q_monthly):
                # try to use 2010 values if they exist
                if np.isnan(q_monthly):
                    q_monthly = monthly_values_2010.get(month, np.nan)
                # otherwise take the average value for each month
                if np.isnan(q_monthly):
                    q_monthly = avg_monthly_values[month]
                # fill missing months with the mean value for the site
                if np.isnan(q_monthly):
                    q_monthly = q_mean
                q.append(q_monthly)
            
            group['q'] = q
            group['q'] = group['q'] * 3785.4  # convert from mgd to cubic m per d

            group['site_no'] = site_no
            group['well_elev_m'] = self.well_elevations_m[site_no]
            group['depth_m'] = self.depths_m[site_no]
            well_botm_depth = self.well_elevations_m[site_no] - self.depths_m[site_no]
            group['x_{0}'.format(self.epsg)] = np.nanmin(group['x_{0}'.format(self.epsg)])
            group['y_{0}'.format(self.epsg)] = np.nanmin(group['y_{0}'.format(self.epsg)])

            # assign a production zone from default dict.  If the bottom of the 
            # well does not fall in a zone, or if the dictionary is empty; then
            # the production zone is assigned 'unnamed'
            production_zone = 'unnamed'
            for prod_name in list(self.prod_zone_top_m):
                prod_zone_top = self.prod_zone_top_m[prod_name][site_no]
                prod_zone_bot = self.prod_zone_bot_m[prod_name][site_no]
                if np.isnan(prod_zone_top) or np.isnan(prod_zone_bot):  # missing zone
                    group['screen_bot_m'] = self.well_elevations_m[site_no] - self.depths_m[site_no]
                    group['screen_top_m'] = self.well_elevations_m[site_no] - self.depths_m[site_no] + self.default_screen_len
                    group['open_int_method'] = 'well depth'
                else:
                    if well_botm_depth < prod_zone_top and well_botm_depth > prod_zone_bot:
                        production_zone = prod_name
                        group['screen_bot_m'] = prod_zone_bot
                        group['screen_top_m'] = prod_zone_top
                        group['open_int_method'] = 'production zone'
                    else:
                        group['screen_bot_m'] = self.well_elevations_m[site_no] - self.depths_m[site_no]
                        group['screen_top_m'] = self.well_elevations_m[site_no] - self.depths_m[site_no] + self.default_screen_len
                        group['open_int_method'] = 'well depth'
            group['production_zone'] = production_zone

            # add aquifer name
            group['aquifer_name'] = self.aquifer_names.get(group["FROM_AQFR_CD"].values[0], 'unnamed')

            cols = ['site_no', 'q', 'q_monthly', 'month', 'well_elev_m', 'depth_m',
                    'screen_bot_m', 'screen_top_m', 'x_{0}'.format(self.epsg), 'y_{0}'.format(self.epsg)]
            all_groups.append(group[cols])

        self.df_swuds = pd.concat(all_groups)
        self.df_swuds['datetime'] = self.df_swuds.index
        self.df_swuds.to_csv(outfile, index=False)
        print('processed SWUDS data written to {0} and in dataframe attribute'.format(outfile))

    @classmethod
    def from_yaml(cls, yamlfile):
        """ Read input and output files from yaml file
        and run all the processing steps in the class to
        produce the processed csv file.

        Parameters
        ----------
        yamlfile: str
            path to a yaml file containing input and output information

        Returns
        -------
        wu: Swuds object
            returns a Swuds object and also generates processed csv file
            specified in the yaml file.

        """
        with open(yamlfile, 'r') as inputs:
            yaml_inputs = yaml.safe_load(inputs)
        
        data_path = yaml_inputs['raw_data']['data_path']
        data_path = Swuds.fix_path(data_path)

        swuds_input = os.path.join(data_path, yaml_inputs['raw_data']['swuds_input'])
        worksheet = yaml_inputs['raw_data']['worksheet']

        outcsv = os.path.join(data_path, yaml_inputs['output']['outcsv'])
        processed_csv = os.path.join(data_path, yaml_inputs['output']['processed_csv'])
        
        dem = os.path.join(data_path,yaml_inputs['rasters']['dem'])
        mc_top = os.path.join(data_path, yaml_inputs['rasters']['mc_top'])
        mc_bot = os.path.join(data_path, yaml_inputs['rasters']['mc_bot'])
        lc_top = os.path.join(data_path, yaml_inputs['rasters']['lc_top'])
        lc_bot = os.path.join(data_path, yaml_inputs['rasters']['lc_bot'])

        meras_shp = os.path.join(data_path, yaml_inputs['shapefiles']['extent'])
        wu_shp = os.path.join(data_path, yaml_inputs['shapefiles']['wells_out'])

        # make a swuds object
        wu = cls(xlsx=swuds_input, sheet=worksheet, csvfile=outcsv)

        # process
        wu.sort_sites()
        wu.reproject()

        wu.apply_footprint(meras_shp, outshp=wu_shp)
        wu.assign_missing_elev(dem)

        for zn in yaml_inputs['zones']:
            zn[1] = Swuds.fix_path(zn[1])
            zn[2] = Swuds.fix_path(zn[2])

        wu.make_production_zones(yaml_inputs['zones'])

        # write results to csv file and return object
        wu.assign_monthly_production(processed_csv)
        return wu
        
    @staticmethod
    def fix_path(data_path):
        """ Convert simple path string with forward slashes
        to a path used by python using os.path.join().  This
        function allows the user to specify a path in the
        yaml file simply, for example:
        d:/home/MAP/source_data/wateruse.csv

        The string is split on '/' and resulting list
        is passed to os.path.join().  If the first entry
        has a colon, then os.path.join(entry[0], os.path.sep)
        is used to specify a windows drive properly.

        Parameters
        ----------
        data_path: str
            path read from yaml file 

        Returns
        -------
        data_path: str
            path built using os.path.join()
        """
        if re.search('/', data_path):
            parts = re.split('/', data_path)
            if re.search(':', parts[0]):
                data_path = os.path.join(parts[0], os.path.sep)
                del parts[0]
                data_path = os.path.join(data_path, *parts)
            else:
                data_path = os.path.join(*parts)
        return(data_path)

    # def overlaps(self, a, b):
    #     """ Return the amount of overlap, in bp
    #     between a and b.

    #     Parameters
    #     ----------
    #     a: 
    #     b: 

    #     Returns
    #     -------
    #     int:
    #         If >0, the number of bp of overlap
    #         If 0,  they are book-ended.
    #         If <0, the distance in bp between them
    #     """

    #     return min(a[1], b[1]) - max(a[0], b[0])


if __name__ == '__main__':

    home = os.getcwd()
    data_path = os.path.join(os.path.dirname(home), 'working_data')
    # swuds_input = os.path.join(data_path, 'LMG-withdrawals-2000-2018.xlsx')
    # worksheet = 'LMG-withdrawals-2000-2018'
    # outcsv = os.path.join(data_path, 'swuds_nonTE.csv')
    # dem = os.path.join(data_path, 'dem_mean_elevs.tif')
    # mc_top = os.path.join(data_path, 'mcaq_surf.tif')
    # mc_bot = os.path.join(data_path, 'mccu_surf.tif')
    # lc_top = os.path.join(data_path, 'lcaq_surf.tif')
    # lc_bot = os.path.join(data_path, 'lccu_surf.tif')

    # meras_shp = os.path.join(data_path, 'MERAS_Extent.shp')
    # wu_shp = os.path.join(data_path, 'WU_points.shp')
    # # make a swuds object
    # # swuds = Swuds(xlsx=swuds_input, sheet=worksheet, csvfile=outcsv)
    # swuds = Swuds(xlsx=None, sheet=None, csvfile=outcsv, cols=None)
    # swuds.sort_sites()
    # swuds.reproject()

    # swuds.apply_footprint(meras_shp, outshp=wu_shp)
    # swuds.assign_missing_elev(dem)

    # mc = ['middle_claiborne', mc_top, mc_bot]
    # lc = ['lower_claiborne', lc_top, lc_bot]
    # swuds.make_production_zones([mc, lc])
    # swuds.assign_monthly_production(os.path.join(data_path, 'processed_swuds.csv'))
    # print(swuds.df_swuds.head())

    yml_file = os.path.join(data_path, 'swuds_input.yml')
    wu = Swuds.from_yaml(yml_file)