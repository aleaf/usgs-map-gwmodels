from collections import defaultdict
from pathlib import Path
import os
import re

import numpy as np
import pandas as pd
from shapely.geometry import Point
import yaml

from gisutils import df2shp, project, raster
from mfsetup.units import convert_length_units, convert_volume_units
from mapgwm.lookups import aq_codes_dict
from mapgwm.utils import cull_data_to_active_area


class Swuds:
    """ Code for preprocessing non-agricultural water use information
    into clean CSV input to MODFLOW setup. Class excludes AQ, IR, and TE
    water use from a swuds excel dataset. Includes logic to fill missing  data,
    as data at many sites are limited to survey years (e.g. 2010 and 2015).
    """

    def __init__(self, xlsx=None, sheet=None, csvfile=None,
                 site_no_col='SITE_NO',
                 x_coord_col='FROM_DEC_LONG_VA',
                 y_coord_col='FROM_DEC_LAT_VA',
                 start_date=None, end_date=None,
                 source_crs=4269, dest_crs=5070,
                 data_length_units='feet',
                 data_volume_units='mgal',
                 model_length_units='meters',
                 default_screen_len=20,
                 cols='default'):
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
        site_no_col : str, optional
            Column name in data with site identifiers,
            by default 'SITE_NO'
        x_coord_col : str, optional
            Column name in data with x-coordinates,
            by default 'x'
        y_coord_col : str, optional
            Column name in data with y-coordinates,
            by default 'y'
        start_date: str
            start time for simulation as string 'yyyy-mm-dd'
        end_date: str
            end time for simulation as string 'yyyy-mm-dd'
        source_crs : obj
            Coordinate reference system of the head observation locations.
            A Python int, dict, str, or :class:`pyproj.crs.CRS` instance
            passed to :meth:`pyproj.crs.CRS.from_user_input`

            Can be any of:
              - PROJ string
              - Dictionary of PROJ parameters
              - PROJ keyword arguments for parameters
              - JSON string with PROJ parameters
              - CRS WKT string
              - An authority string [i.e. 'epsg:4326']
              - An EPSG integer code [i.e. 4326]
              - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
              - An object with a `to_wkt` method.
              - A :class:`pyproj.crs.CRS` class

            By default, epsg:4269
        dest_crs : obj
            Coordinate reference system of the model. Same input types
            as ``source_crs``.
            By default, epsg:5070
        data_length_units : str; 'meters', 'feet', etc.
            Units of lengths in data (elevations, depths, etc.)
            by default, 'feet'
        data_volume_units : str; 'mgd', 'ft3', etc
            Volumetric unit of pumping rates
            by default, 'mgd' (million gallons per day)
        model_length_units : str; 'meters', 'feet', etc.
            Length units of model.
            by default, 'meters'

        Attributes
        ----------
        df: pandas dataframe
            pandas dataframe read from spreadsheet or csv.  Manipulated by other
            methods of the class
        aquifer_names: dict
            dictionary of aquifer names keyed by NWIS codes, read using import 
            statement from mapgwm.lookups
        regional_aquifers: dict
            dictionary of regional aquifers keyed by NWIS codes, read using import
            statement from mapgwm.lookups
        crs: int
            pyproj.crs.CRS instance describing the output coordinate reference system
        monthly_cols: list of str
            list of column names for monthly values
        start_date: str
            start time for simulation as string 'yyyy-mm-dd'
        end_date: str
            end time for simulation as string 'yyyy-mm-dd'
        default_screen_len: float
            default screen length in meters
        locations: dict
            dictionary of x,y locations keyed by SITE_NO, added in reproject method
        depths: dict
            dictionary of depth keyed by SITE_NO
        well_elevations: dict
            dictionary of well elevations keyed by SITE_NO
        prod_zone_top: defaultdict(dict)
            defaultdict with production zone top (in meters) for each well
            First key is the production zone name, and second is the SITE_NO
            For example  self.prod_zone_top['lower_claiborne']['WEL001'] = top_elev
        prod_zone_bot: defaultdict(dict)
            defaultdict with production zone bottom (in meters) for each well
            First key is the production zone name, and second is the SITE_NO
            For example  self.prod_zone_bot['lower_claiborne']['WEL001'] = bot_elev
        
        """

        # set some attributes
        self.monthly_cols = ['JAN_VAL', 'FEB_VAL', 'MAR_VAL', 'APR_VAL', 'MAY_VAL',
                             'JUN_VAL', 'JUL_VAL', 'AUG_VAL', 'SEP_VAL', 'OCT_VAL',
                             'NOV_VAL', 'DEC_VAL']
        self.aquifer_names = aq_codes_dict['aquifer_code_names']
        self.regional_aquifers = aq_codes_dict['regional_aquifer']
        self.start_date = start_date
        self.end_date = end_date
        self.default_screen_len = default_screen_len
        self.source_crs = source_crs
        self.dest_crs = dest_crs
        self.data_length_units = data_length_units
        self.data_volume_units = data_volume_units
        self.model_length_units = model_length_units
        self.prod_zone_top = defaultdict(dict)
        self.prod_zone_bot = defaultdict(dict)
        self.locations = dict()
        self.site_no_col = site_no_col

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
            raise ValueError('both xlsx and csvfile cannont be None for SWUDS object')
        elif xlsx is None:
            self.df = pd.read_csv(csvfile, usecols=usecols)
        else:
            self.df = pd.read_excel(xlsx, sheet_name=sheet, usecols=usecols)
            if csvfile is not None:
                self.df.to_csv(csvfile, index=False)

        self.df.columns = self.df.columns.str.upper()
        # remove trailing spaces from code names
        if 'FROM_AQFR_CD' in self.df.columns:
            self.df["FROM_AQFR_CD"] = self.df["FROM_AQFR_CD"].str.strip()

        # make depths floats
        if 'FROM_WELL_DEPTH_VA' in self.df.columns:
            self.df['SCREEN_BOT'] = pd.to_numeric(self.df['FROM_WELL_DEPTH_VA'], errors='coerce')
        else:
            self.df['SCREEN_BOT'] = np.nan
        if 'FROM_ALT_VA' in self.df.columns:
            self.df['FROM_ALT_VA'] = pd.to_numeric(self.df['FROM_ALT_VA'], errors='coerce')

        # make dictionaries
        length_conversion = convert_length_units(data_length_units, model_length_units)
        self.depths = dict(list(zip(self.df['SITE_NO'], self.df['SCREEN_BOT'] * length_conversion)))
        self.well_elevations = dict(zip(self.df['SITE_NO'], self.df['FROM_ALT_VA'] * length_conversion))

        self.sort_sites(primarysort=site_no_col)
        # Best to reproject on init so that we know the points are in the dest_crs
        self.reproject(x_coord_col=x_coord_col, y_coord_col=y_coord_col, key=site_no_col)

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
            self.df = self.df.sort_values([primarysort, secondarysort]).groupby([primarysort]).last().reset_index()
        else:
            self.df = self.df.sort_values([primarysort]).groupby([primarysort]).last().reset_index()

    def reproject(self, x_coord_col='FROM_DEC_LONG_VA', y_coord_col='FROM_DEC_LAT_VA', key='SITE_NO'):
        """ Reproject from self.source_crs to self.dest_crs using gisutils

        Parameters
        ----------
        x_coord_col : str, optional
            Column name in data with x-coordinates,
            by default 'x'
        y_coord_col : str, optional
            Column name in data with y-coordinates,
            by default 'y'
        key: str
            key for the dictionary made, defaults to SITE_NO
        """
        x_reprj, y_reprj = project(zip(self.df[x_coord_col], self.df[y_coord_col]),
                                   self.source_crs, self.dest_crs)
        self.df['x'] = x_reprj
        self.df['y'] = y_reprj
        self.df['geometry'] = [Point(x, y) for x, y in zip(x_reprj, y_reprj)]

        # drop entries if no location information
        self.df.dropna(subset=['x', 'y'], axis=0, inplace=True)

        # make dictionary of location 
        self.locations = dict(list(zip(self.df[key], list(zip(self.df['x'], self.df['y'])))))

    def apply_footprint(self, active_area, active_area_id_column=None,
                     active_area_feature_id=None,):
        """Keep sites in the df pandas dataframe that fall
        into the passed bounding shapefile polygon. Requires
        that df dataframe has a Point geometry column as assigned
        in the reproject method.

        Parameters
        ----------
        active_area: str
            path to shapefile with footprint for current analysis
        active_area_id_column : str, optional
            Column in active_area with feature ids.
            By default, None, in which case all features are used.
        active_area_feature_id : str, optional
            ID of feature to use for active area
            By default, None, in which case all features are used.
        outshp: str
            optional path to output shapefile with points within the footprint
        """
        self.df = cull_data_to_active_area(self.df, active_area,
                                           active_area_id_column,
                                           active_area_feature_id,
                                           data_crs=self.dest_crs)

    def assign_missing_elevs(self, top_raster, dem_units='meters'
                             ):
        """ Use the top of model raster, or land-surface raster,
        to assign the elevation for points where elevation is missing.
        
        Parameters
        ----------
        top_raster: str
            path to raster data set with land surface or model top elevation, used
            to assign missing values to water-use points
        elev_field: str
            field in df with elevation data, default is 'FROM_ALT_VA'
        """

        no_elev = self.df['FROM_ALT_VA'].isnull()
        x_no_elev = self.df.loc[no_elev, 'x'].values
        y_no_elev = self.df.loc[no_elev, 'y'].values
        elevs = raster.get_values_at_points(top_raster,
                                            x=x_no_elev,
                                            y=y_no_elev,
                                            points_crs=self.dest_crs
                                            )
        elevs *= convert_length_units(dem_units, self.model_length_units)
        self.df.loc[no_elev, 'FROM_ALT_VA'] = elevs
        assert not self.df['FROM_ALT_VA'].isnull().any()
        self.well_elevations = dict(zip(self.df['SITE_NO'], self.df['FROM_ALT_VA']))

    def make_production_zones(self, production_zones, default_elevation_units='feet'):
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
        #if isinstance(zonelist[0], str):
        #    zonelist = [zonelist]
        key = self.site_no_col

        # get tops and bottoms of estimated production intervals at each well
        # make dictionaries to lookup by well
        for name, info in production_zones.items():
            top_raster, botm_raster, *units = info
            units = units[0] if units else default_elevation_units
            x = self.df['x'].values
            y = self.df['y'].values
            length_unit_conversion = convert_length_units(units, self.model_length_units)
            top_elevations = raster.get_values_at_points(top_raster,
                                                 x=x, y=y) * length_unit_conversion
            self.prod_zone_top[name] = dict(zip(self.df[key], top_elevations))
            botm_elevations = raster.get_values_at_points(botm_raster,
                                                          x=x, y=y) * length_unit_conversion
            self.prod_zone_bot[name] = dict(zip(self.df[key], botm_elevations))
            self.df['{}_top'.format(name)] = top_elevations
            self.df['{}_botm'.format(name)] = botm_elevations

    def assign_monthly_production(self, outfile='processed_swuds.csv'):
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
            idx = self.df.loc[self.df[c].isnull()].index.values
            self.df.loc[idx, c] = self.df.loc[idx, 'ANNUAL_VAL']

        # pull out groundwater sites that are not IR, AQ or TE
        self.df = self.df.loc[
            (self.df['WATER_CD'] == 'GW') &
            ~(self.df['FROM_NAT_WATER_USE_CD'] == 'IR') &
            ~(self.df['FROM_NAT_WATER_USE_CD'] == 'AQ') &
            ~(self.df['FROM_NAT_WATER_USE_CD'] == 'TE')
            ]

        # reshape dataframe to have monthly values in same column
        stacked = pd.DataFrame(self.df[self.monthly_cols].stack())
        stacked.reset_index(inplace=True)
        stacked.rename(columns={'level_1': 'month',
                                0: 'q_monthly'}, inplace=True)
        stacked.q_monthly = stacked.q_monthly
        stacked.index = stacked.level_0
        stacked = stacked.join(self.df)
        keep_cols = [c for c in stacked.columns if c not in self.monthly_cols]
        stacked = stacked[keep_cols]
        month = {name: i + 1 for i, name in enumerate(self.monthly_cols)}
        dates = ['{}-{:02d}'.format(year, month[month_column_name])
                for year, month_column_name in zip(stacked.YEAR, stacked.month)]
        stacked['datetime'] = pd.to_datetime(dates)
        stacked.sort_values(by=['SITE_NO', 'datetime'], inplace=True)

        # set start and end dates if not already set
        if self.start_date is None:
            self.start_date = stacked.datetime.min()
        if self.end_date is None:
            self.end_date = stacked.datetime.max()

        groups = stacked.groupby('SITE_NO')
        all_groups = []
        for site_no, group in groups:
            group = group.copy()
            group.index = pd.to_datetime(group['datetime'])
            start_date = pd.Timestamp(self.start_date)
            end_date = pd.Timestamp(self.end_date)

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
            #group['q'] = group['q'] * 3785.4  # convert from mgd to cubic m per d
            group['q'] = group['q'] * convert_volume_units(self.data_volume_units,
                                                           self.model_length_units)

            group['site_no'] = site_no
            group['well_elev'] = self.well_elevations[site_no]
            group['depth'] = self.depths[site_no]
            well_botm_depth = self.well_elevations[site_no] - self.depths[site_no]
            group['x'] = np.nanmin(group['x'])
            group['y'] = np.nanmin(group['y'])

            # assign a production zone from default dict.  If the bottom of the 
            # well does not fall in a zone, or if the dictionary is empty; then
            # the production zone is assigned 'unnamed'
            production_zone = 'unnamed'
            for prod_name in self.prod_zone_top.keys():
                prod_zone_top = self.prod_zone_top[prod_name][site_no]
                prod_zone_bot = self.prod_zone_bot[prod_name][site_no]
                if np.isnan(prod_zone_top) or np.isnan(prod_zone_bot):  # missing zone
                    group['screen_bot'] = self.well_elevations[site_no] - self.depths[site_no]
                    group['screen_top'] = self.well_elevations[site_no] - self.depths[site_no] + self.default_screen_len
                    group['open_int_method'] = 'well depth'
                else:
                    if well_botm_depth < prod_zone_top and well_botm_depth > prod_zone_bot:
                        production_zone = prod_name
                        group['screen_bot'] = prod_zone_bot
                        group['screen_top'] = prod_zone_top
                        group['open_int_method'] = 'production zone'
                    else:
                        group['screen_bot'] = self.well_elevations[site_no] - self.depths[site_no]
                        group['screen_top'] = self.well_elevations[site_no] - self.depths[site_no] + self.default_screen_len
                        group['open_int_method'] = 'well depth'
            group['production_zone'] = production_zone

            # add aquifer name
            group['aquifer_name'] = self.aquifer_names.get(group["FROM_AQFR_CD"].values[0], 'unnamed')

            cols = ['site_no', 'q', 'q_monthly', 'month', 'well_elev', 'depth',
                    'screen_bot', 'screen_top', 'x', 'y']
            all_groups.append(group[cols])

        self.df = pd.concat(all_groups)
        self.df['start_datetime'] = self.df.index  # start date of each pumping period
        if outfile is not None:
            outfile = Path(outfile)
            self.df.to_csv(outfile, index=False)
            print('processed SWUDS data written to {0} and in dataframe attribute'.format(outfile))
            self.df['geometry'] = [Point(x, y) for x, y in zip(self.df.x,
                                                               self.df.y)]
            # write only unique pumping values to shapefile
            to_shapefile = self.df.groupby(['site_no', 'q']).first().reset_index()
            shapefile = outfile.with_suffix('.shp')
            df2shp(to_shapefile, shapefile, crs=self.dest_crs)

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
        wu.assign_missing_elevs(dem)

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


def preprocess_swuds(swuds_input, worksheet, csv_input=None,
                     dem=None, dem_units='meters',
                     start_date=None,
                     end_date=None,
                     active_area=None,
                     active_area_id_column=None,
                     active_area_feature_id=None,
                     site_no_col='SITE_NO',
                     x_coord_col='FROM_DEC_LONG_VA',
                     y_coord_col='FROM_DEC_LAT_VA',
                     production_zones=None,
                     estimated_production_surface_units='feet',
                     source_crs=4269,
                     dest_crs=5070,
                     data_length_units='feet',
                     data_volume_units='mgal',
                     model_length_units='meters',
                     outfile=None):
    """Preprocess water use data from the USGS Site-Specific Water Use Database (SWUDS).

    * reproject data to a destination CRS `dest_crs`)
    * cull data to an area of interest (`active_area`)
    * assign any missing wellhead elevations from a DEM
    * if input data do not have information on the well screen intervals;
      sample screen tops and bottoms from raster surfaces bounding
      an estimated production zone (e.g. `estimated_production_zone_top`). Well
      bottom information is used to discriminate between multiple production zones.
    * reindex the data to continous monthly values extending from `start_date`
      to `end_date`. Typically, these would bracket the time period for which
      the pumping should be simulated in a model. For example, the earliest data
      may be from 2010, but if the model starts in 2008, it may be appropriate to
      begin using the 2010 rates then (``start_date='2008'``). If no start or end
      date are given, the first and last years of pumping in `data` are used.
    * fill empty months using 2010 data (the most complete survey year) if available,
      otherwise use the average value for the site.
    * backfill any remaining empty months going back to the `start_date` in the same way
    * write processed data to a CSV file and shapefile of the same name

    Parameters
    ----------
    swuds_input: str
        Excel spreadsheet of SWUDs data, in a format readable by
        :func:`pandas.read_excel`. If xlsx file is passed, then a
        selected worksheet from it will be converted to a csv file unless
        the csvfile parameter is specified as None.
    worksheet : str
        Worksheet in `swuds_input` to read.
    csvfile: str, optional
        Path to csv file with data (if xlsx if None) or to a csvfile
        that is created from the selected worksheet.  If xlsx is None,
        then csvfile must be provided.
    sheet: str
        Name of worksheet in xlsx to be read, ignored if xlsx is None.
    dem : str, optional
        DEM raster of the land surface. Used for estimating missing wellhead elevations.
        Any reprojection to dest_crs is handled automatically, assuming
        the DEM raster has CRS information embedded (arc-ascii grids do not!)
        By default, None.
    dem_units : str, {'feet', 'meters', ..}
        Units of DEM elevations, by default, 'meters'
    start_date : str
        Start date for pumping rates. If earlier than the dates in `data`,
        pumping rates will be backfilled to this date.
    end_date : str
        End date for pumping rates. If later than the dates in `data`,
        pumping rates will be forward filled to this date.
    active_area : str
        Shapefile with polygon to cull observations to. Automatically reprojected
        to dest_crs if the shapefile includes a .prj file.
        by default, None.
    active_area_id_column : str, optional
        Column in active_area with feature ids.
        By default, None, in which case all features are used.
    active_area_feature_id : str, optional
        ID of feature to use for active area
        By default, None, in which case all features are used.
    site_no_col : str, optional
        Column name in data with site identifiers,
        by default 'SITE_NO'
    x_coord_col : str, optional
        Column name in data with x-coordinates,
        by default 'FROM_DEC_LONG_VA'
    y_coord_col : str, optional
        Column name in data with y-coordinates,
        by default 'FROM_DEC_LAT_VA'
    production_zones : dict
        Dictionary of production zone tops and bottoms, and optionally,
        production zone surface elevation units, keyed by an abbreviated name.

        Example::

            production_zones={'mrva': (test_data_path / 'swuds/rasters/pz_MCAQP_top.tif',
                                       test_data_path / 'swuds/rasters/pz_MCAQP_bot.tif',
                                    ),
                              'mcaq': (test_data_path / 'swuds/rasters/pz_MCAQP_top.tif',
                                       test_data_path / 'swuds/rasters/pz_MCAQP_bot.tif',
                                       'feet')
                              }

        Where zone 'mrva' has a top and bottom raster, but no units assigned,
        and zone 'mcaq' has a top and bottom raster and units. For zones
        with no units, the `estimated_production_surface_units` will be used.

    estimated_production_surface_units : str, {'meters', 'ft', etc.}
        Length units of elevations in estimated production surface rasters.
        by default, 'feet'
    source_crs : obj
        Coordinate reference system of the head observation locations.
        A Python int, dict, str, or :class:`pyproj.crs.CRS` instance
        passed to :meth:`pyproj.crs.CRS.from_user_input`

        Can be any of:
          - PROJ string
          - Dictionary of PROJ parameters
          - PROJ keyword arguments for parameters
          - JSON string with PROJ parameters
          - CRS WKT string
          - An authority string [i.e. 'epsg:4326']
          - An EPSG integer code [i.e. 4326]
          - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
          - An object with a `to_wkt` method.
          - A :class:`pyproj.crs.CRS` class

        By default, epsg:4269

    dest_crs : obj
        Coordinate reference system of the model. Same input types
        as ``source_crs``.
        By default, epsg:5070
    data_length_units : str; 'meters', 'feet', etc.
        Units of lengths in data (elevations, etc.)
        by default, 'meters'
    data_volume_units : str; 'mgd', 'ft3', etc
        Volumetric unit of pumping rates
        by default, 'mgd' (million gallons per day)
    model_length_units : str; 'meters', 'feet', etc.
        Length units of model.
        by default, 'meters'
    outfile : str
        Path for output file. A shapefile of the same name is also written.
        If None, no output file is written. By default, None

    Returns
    -------
    wu : :class:`mapgwm.swuds.Swuds` instance

    Notes
    -----
    * time units for SWUDs data and model are assumed to be days

    """
    # make a swuds object and process it
    wu = Swuds(xlsx=swuds_input, sheet=worksheet, csvfile=csv_input,
               start_date=start_date,
               end_date=end_date,
               site_no_col=site_no_col,
               x_coord_col=x_coord_col,
               y_coord_col=y_coord_col,
               source_crs=source_crs,
               dest_crs=dest_crs,
               data_length_units=data_length_units,
               data_volume_units=data_volume_units,
               model_length_units=model_length_units,
               )

    # cull the data to the model area, if provided
    if active_area is not None:
        wu.apply_footprint(active_area=active_area,
                           active_area_id_column=active_area_id_column,
                           active_area_feature_id=active_area_feature_id)
    wu.assign_missing_elevs(dem, dem_units=dem_units
                            )

    #mc = ['middle_claiborne', mc_top, mc_bot]
    #lc = ['lower_claiborne', lc_top, lc_bot]
    if production_zones is not None:
        wu.make_production_zones(production_zones,
                                 default_elevation_units=estimated_production_surface_units)
    wu.assign_monthly_production(outfile=outfile)
    return wu


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
    # swuds.assign_missing_elevs(dem)

    # mc = ['middle_claiborne', mc_top, mc_bot]
    # lc = ['lower_claiborne', lc_top, lc_bot]
    # swuds.make_production_zones([mc, lc])
    # swuds.assign_monthly_production(os.path.join(data_path, 'processed_swuds.csv'))
    # print(swuds.df.head())

    yml_file = os.path.join(data_path, 'swuds_input.yml')
    wu = Swuds.from_yaml(yml_file)