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
        # merge the individual dataframes if there is more than one
        if len(tes) > 1:
            self.df_swuds = pd.concat(tes, axis=1)
            # get rid of repeated columns from concat, duplicated() returns true
            # for repeated names, so select the columns not-duplicated()
            self.df_swuds = self.df_swuds.loc[:, ~self.df_swuds.columns.duplicated()]  
        else:
            self.df_swuds = tes[0]


    def assign_monthly_production(self, outfile, year_cols):
        """ Assign production wells for TE water use.  Overrides the method
        from Swuds() class.  The TE model does not include well information, 
        so the production zones assigned to the object are used with the
        default screen length to define the pumping interval.  If a location
        does not have have production zone information, then nan is assigned
        and the user will have to assign tops and bottoms to the dataframe
        and final csv.  

        Production is given in cubic m per day.

        todo:  add unit conversion parameter so other units can be used?

        Parameters
        ----------
        outfile: str
            path to final processed monthly TE water-use file.
        year_cols: list of str
            list with columns to be gathered to assign monthly values
            Usually these are from 5-year TE model runs.  Should contain
            a 4-digit year data that can be matched in each string. These
            will correpond to the names assigned in the constructor.
        """

        # pull out groundwater sites 
        self.df_swuds = self.df_swuds.loc[self.df_swuds['src_code'] == 'GW']

        # get years from year_cols
        years = []
        col_dict = dict()
        for col in year_cols:
            m = re.search('(\d\d\d\d)', col)
            if m:
                years.append(int(m[1]))
                col_dict[int(m[1])] = col

        firstyear = pd.to_datetime('{0}-1-1'.format(min(years)))
        midpoints = []
        for i in range(0, len(years)-1):
            y, r = divmod((years[i]+years[i+1]), 2)
            days = r/2. * 365.
            midpoints.append(pd.to_datetime('{0}-1-1'.format(y)) + pd.to_timedelta(days, unit='days'))

        df = self.df_swuds.copy()
        
        df['start_date'] = firstyear
        
        groups = df.groupby('plant_code')
        all_groups = []
        for plant_code, group in groups:
            group = group.copy()
            group.index = group['start_date']
            start_date = pd.Timestamp(self.sim_start_dt)
            end_date = pd.Timestamp(self.sim_end_dt)
            parts = list(map(int, re.split('-', self.sim_start_dt)))
            start_datetime = pd.to_datetime(self.sim_start_dt)

            monthly_values=dict()
            for y in years:
                monthly_values[y] = group.loc[group.start_date == firstyear,col_dict[y]].values[0]

            all_dates = pd.date_range(start_date, end_date, freq='MS')
            group = group.reindex(all_dates)

            # fill empty dates
            mask = ((group.index >= start_datetime) & (group.index < midpoints[0]))
            group.loc[mask, 'q_m3d'] = monthly_values[years[0]] * 3785.4   # convert to m3 per day

            for i in range(1, len(midpoints)):
                mask = (group.index >= datetime.date(midpoints[i-1]) & (group.index < midpoints[i]))
                group.loc[mask, 'q_m3d'] = monthly_values[years[i]] * 3785.4

            # take the first key in the production zones if more than one
            prod_zone = list(self.prod_zone_top_m.keys())[0]
            
            group['screen_top_m'] = self.prod_zone_top_m[prod_zone][plant_code]
            group['screen_bot_m'] = self.prod_zone_top_m[prod_zone][plant_code]
            group['x_5070'] = np.nanmin(group['x_5070'])
            group['y_5070'] = np.nanmin(group['y_5070'])
            group['plant_code'] = plant_code

            cols = ['plant_code', 'q_m3d', 'screen_bot_m',
                    'screen_top_m', 'x_5070', 'y_5070']
            all_groups.append(group[cols])

        self.df_swuds = pd.concat(all_groups)
        self.df_swuds['datetime'] = self.df_swuds.index
        self.df_swuds.to_csv(outfile, index=False)
        print('processed TE SWUDS data written to {0} and in dataframe attribute'.format(outfile))

    @classmethod
    def from_yaml(cls, yamlfile):
        """ Read input and output files from yaml file
        and run all the processing steps in the class to
        produce the processed csv file.  Overrided method
        from Swuds class to do the steps for 
        ThermalUse object. 

        Parameters
        ----------
        yamlfile: str
            path to a yaml file containing input and output information

        Returns
        -------
        te: ThermalUse object
            returns a ThermalUse object and also generates processed csv file
            specified in the yaml file.

        """
        with open(yamlfile, 'r') as inputs:
            yaml_inputs = yaml.safe_load(inputs)
        
        data_path = ThermalUse.fix_path(yaml_inputs['raw_data']['data_path'])

        te_input = []
        for f in yaml_inputs['raw_data']['te_input']:
            te_input.append(os.path.join(data_path, f))
        worksheets = yaml_inputs['raw_data']['worksheet']

        outcsvs = []
        for f in yaml_inputs['output']['outcsv']:
            outcsvs.append(os.path.join(data_path, f))
        
        processed_csv = os.path.join(data_path, yaml_inputs['output']['processed_csv'])
        
        dem = os.path.join(data_path,yaml_inputs['rasters']['dem'])
        mc_top = os.path.join(data_path, yaml_inputs['rasters']['mc_top'])
        mc_bot = os.path.join(data_path, yaml_inputs['rasters']['mc_bot'])
        lc_top = os.path.join(data_path, yaml_inputs['rasters']['lc_top'])
        lc_bot = os.path.join(data_path, yaml_inputs['rasters']['lc_bot'])

        meras_shp = os.path.join(data_path, yaml_inputs['shapefiles']['extent'])
        te_shp = os.path.join(data_path, yaml_inputs['shapefiles']['wells_out'])

        # make a te object
        te = cls(xlsx=te_input, 
                    sheet=worksheets, 
                    csvfile=outcsvs, 
                    cols=yaml_inputs['raw_data']['cols'], 
                    skiprows=yaml_inputs['raw_data']['skiprows'], 
                    newheaders=yaml_inputs['raw_data']['newheaders'], 
                    epsg=5070)

        # process
        te.sort_sites(primarysort=yaml_inputs['options']['primarysort'])
        te.reproject(long=yaml_inputs['options']['long'], 
                     lat=yaml_inputs['options']['lat'],
                     key=yaml_inputs['options']['primarysort'])

        te.apply_footprint(meras_shp, outshp=te_shp)

        for zn in yaml_inputs['zones']:
            zn[1] = Swuds.fix_path(zn[1])
            zn[2] = Swuds.fix_path(zn[2])

        te.make_production_zones(yaml_inputs['zones'],
                                 key=yaml_inputs['options']['primarysort'])

        # write results to csv file and return object
        te.assign_monthly_production(processed_csv, 
                                     year_cols=yaml_inputs['options']['year_cols'])
        return te

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

    # write the combined dataframe 
    te.df_swuds.to_csv(os.path.join(data_path, 'TE_combined.csv'))

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
    
    mc = ['middle_claiborne', mc_top, mc_bot]
    te.make_production_zones([mc], key='plant_code')
    
    te.assign_monthly_production(os.path.join(data_path, 'processed_thermal.csv'), 
                                 year_cols=['2010_q_mgd', '2015_q_mgd'])
    print(te.df_swuds.head())

    # try the yaml method
    yml_file = os.path.join(data_path, 'thermal_input.yml')
    te = ThermalUse.from_yaml(yml_file)
    print(te.df_swuds.head())