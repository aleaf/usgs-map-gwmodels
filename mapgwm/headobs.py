import os
import numpy as np
import pandas as pd
from shapely.geometry import Point
import matplotlib.pyplot as plt
from gisutils import shp2df, df2shp, project
from mfsetup.obs import make_obsname
from mfsetup.units import convert_length_units


def makedirs(path):
    dirs = [os.path.join(path, 'shps/'),
            os.path.join(path, 'figures/')]
    for folder in dirs:
        if not os.path.isdir(folder):
            os.makedirs(folder)


def read_metadata(files, usecols=None):
    """

    :param zmetapath:file path for metadata files
    :return: indf : pandas dataframe with metadata
    """
    dfs = []
    for f in files:
        df = pd.read_csv(f, usecols=usecols,
                         sep='\t', comment='#', low_memory=False)
        dfs.append(df)
    df = pd.concat(dfs, sort=True)
    return df


def get_layer(botm, i, j, elev, nlay):
    """Return the layers for elevations at i, j locations.

    Parameters
    ----------
    botm : bottom array from modflow model
    i : scaler or sequence
        row index (zero-based)
    j : scaler or sequence
        column index
    elev : scaler or sequence
        elevation (in same units as model)
    nlay : number of layers
    Returns
    -------
    k : np.ndarray (1-D) or scalar
        zero-based layer index
    """
    def to_array(arg):
        if not isinstance(arg, np.ndarray):
            return np.array([arg])
        else:
            return arg

    i = to_array(i)
    j = to_array(j)
    elev = to_array(elev)
    botms = botm[:, i, j].tolist()
    layers = np.sum(((botms - elev) > 0), axis=0)
    # force elevations below model bottom into bottom layer
    layers[layers > nlay - 1] = nlay - 1
    layers = np.atleast_1d(np.squeeze(layers))
    if len(layers) == 1:
        layers = layers[0]
    return layers


def get_active_area(shapefile, name_col='desc', 
                    buffer=10000.):
    """Read a polygon from a shapefile and add a buffer.

    Parameters
    ----------
    shapefile : [type]
        [description]
    name_col : str, optional
        [description], by default 'desc'
    buffer : [type], optional
        [description], by default 10000.

    Returns
    -------
    [type]
        [description]
    """    
    df = shp2df(shapefile)
    #df.index = df[name_col]
    # take the first polygon by default
    return df.geometry.values[0].buffer(buffer)
    
    
def get_header_length(sitefile, col0='SITE_BADGE'):
    
    with open(sitefile) as src:
        for i, line in enumerate(src):
            if '#' not in str(line) and col0 in str(line):
                return i
            
            
def get_data(data_file, metadata_file):
    
    # inputs
    # rename columns for readability, shapefile 10 char. limit
    # and to conform with mfsetup naming conventions
    # for now, read all of the columns in data_file
    data_usecols = None
    data_column_renames = {# columns in 4/14/20 "monthly rollup"
                           'SITE_BADGE': 'site_no',
                           'DEC_LAT_VA': 'lat',
                           'DEC_LONG_VA': 'lon',
                           'MONTHLY_LEV_ALT_VA': 'head_m',
                           'MONTHLY_STDEV_LEV_ALT_VA': 'head_std_m',
                           'MONTHLY_LAST_LEV_ALT_VA': 'last_head_m',
                           'MONTHLY_LEV_COUNT': 'n',
                           # columns in previous (statistical model output) version
                           # (in case they're needed again)
                           'LEV_DT': 'datetime',
                           'EST_LEV_ALT_STDEV': 'head_std_m',
                           'ALT_VA': 'well_el_m',
    }
    # read only these columns from metadata_file
    metadata_usecols = [
                        'SITE_BADGE',
                        'WELL_DEPTH_VA',
                        'OPEN_TOP_VA',
                        'OPEN_BOTTOM_VA',
                        'AQFR_CD',
                        'NAT_AQFR_CD',
                        'ALT_VA'
                        ]
    # rename these columns
    metadata_column_renames = {'SITE_BADGE': 'site_no',
                               'WELL_DEPTH_VA': 'well_depth',
                               'OPEN_TOP_VA': 'screen_top',
                               'OPEN_BOTTOM_VA': 'screen_botm',
                               'ALT_VA': 'well_el_m'
    }
    
    # names associated with aquifer codes
    # from https://help.waterdata.usgs.gov/aqfr_cd
    aqfr_cd_name = {'124SPRT': 'Sparta Sand', 
                    '110ALVM': 'Quaternary Alluvium', 
                    '124MUWX': 'Meridian-Upper Wilcox Aquifer', 
                    '111ALVM': 'Holocene Alluvium', 
                    '124MMPS': 'Memphis Sand', 
                    '112TRRC': 'Terrace Deposits',
                    '211MCNR': 'McNairy Formation', 
                    '124WLCXG': 'Wilcox Group', 
                    '12405MP': '500-foot Sand (Memphis Aquifer)', 
                    '112MRVA': 'MRVA', 
                    '112MRVAA': 'MRVA',
                    '112VLTR': 'Valley Trains, Pleistocene',
                    '112ACFL': 'Atchafalaya Aquifer',
                    '112ORVA': 'Ouachita River Alluvial Aquifer',
                    '110QRNR': 'Quaternary System',
                    '112ABVL': 'Abbeville Aquifer, Pleistocene Age',
                    '112PRIR': 'Prairie Formation',
                    '112LOSS': 'Loess',
                    '111HLCN': 'Holocene Series',
                    '112ALVL': 'Alluvial Aquifers, Undifferentiated',
                    '124WLCX': 'Wilcox Group', 
                    '124FRPL': 'Fort Pillow Sand of Wilcox Group',
                    'N100MSRVVL': 'MRVA'
    }
    
    # read in the data
    data_skiprows = get_header_length(data_file)
    data = pd.read_csv(data_file, usecols=data_usecols, sep='\t', 
                       skiprows=data_skiprows)
    data.rename(columns=data_column_renames, inplace=True)
    data.columns = data.columns.str.lower()
    if 'datetime' not in data.columns:
        datetimes = ['{}-{:02d}'.format(year, month) 
                     for year, month in zip(data.year, data.month)]
        data['datetime'] = pd.to_datetime(datetimes)
    


    # read in the metadata
    metadata_skiprows = get_header_length(metadata_file)
    metadata = pd.read_csv(metadata_file, usecols=metadata_usecols, sep='\t', 
                           skiprows=metadata_skiprows)
    metadata.rename(columns=metadata_column_renames, inplace=True)
    metadata.set_index('site_no', inplace=True)
    metadata.columns = metadata.columns.str.lower()
    
    well_el_m = dict(zip(metadata.index, metadata['well_el_m']))
    data['well_el_m'] = [well_el_m[site_no] for site_no in data.site_no]

    print('metatadata has no column called "aquifer, so must create it from local or national aquifer codes"')
    #create aquifer column in metadata
    metadata.fillna(value={'aqfr_cd':'unspecified', 'nat_aqfr_cd':'unspecified'}, inplace=True)
    metadata.loc[metadata.aqfr_cd != 'unspecified', 'aquifer'] = metadata.aqfr_cd
    metadata.loc[(metadata.aqfr_cd == 'unspecified') & (metadata.nat_aqfr_cd != 'unspecified'), 'aquifer'] = metadata.nat_aqfr_cd
    metadata.loc[(metadata.aqfr_cd == 'unspecified') & (metadata.nat_aqfr_cd == 'unspecified'), 'aquifer'] = 'unspecified' 
    
    #set(data.aquifer).difference(set(aqfr_cd_name.keys()))
    metadata['aquifer_name'] = [aqfr_cd_name[aqfrcdi] for aqfrcdi in metadata.aquifer]
    return data, metadata


def preprocess_headobs(data, metadata, start_date='1998-04-01',
                       src_epsg=4269, dest_epsg=5070,
                       data_length_units='meters',
                       active_area=None, active_area_epsg=5070,
                       aoi=None,
                       max_obsname_len=13,
                       outfile='../source_data/observations/head_obs/preprocessed_head_obs.csv'):
    """

    Parameters
    ----------
    data : DataFrame
        Head observation data
    metadata : DataFrame
        Well information
    start_date : str (YYYY-mm-dd)
        Simulation start date (cull observations before this date)
    src_epsg : int
        EPSG code defining projection for source data
    dest_epsg : int
        EPSG code defining projection of model
    data_length_units : str; 'meters', 'feet', etc.
        Length units of head observations.
    active_area : polygon
    active_area_epsg : int
        EPSG code defining CRS for active area polygon
    aoi : dict  (name: shapefile path)
        Areas of interest. Dictionary of shapefile paths
        keyed by area names as they will in observation group names
        (e.g. 'mscha' for MS composite hydrograph area').
    max_obsname_len : int
        Maximum length for observation name prefix. Default of 13
        allows for a PEST obsnme of 20 characters or less with
        <prefix>_yyyydd or <prefix>_<per>d<per>
        (e.g. <prefix>_2d1 for a difference between stress periods 2 and 1)
    outfile : str
        Where output file will be written

    Returns
    -------

    """

    df = data.copy()
    # multiplier to convert input length units to model units
    unit_conversion = convert_length_units(data_length_units, 'meters')

    # outputs
    outpath, filename = os.path.split(outfile)
    makedirs(outpath)
    outname, ext = os.path.splitext(outfile)
    out_info_csvfile = outname + '_info.csv'
    out_data_csvfile = outfile
    out_plot = os.path.join(outpath, 'open_interval_lengths.pdf')
    out_shapefile = outname + '_info.shp'

    # set the starting and ending dates here
    stdate = pd.Timestamp(start_date)

    # convert to datetime; drop the timestamps
    df['datetime'] = pd.to_datetime(df.datetime).dt.normalize()

    # trim to the time range
    df = df.loc[(df.datetime >= stdate)]

    # make two tables- one with well info and other with just data, site numbers, times and aquifer
    #aquifer = dict(zip(df.site_no, df.aquifer))
    
    # well info
    # collapse dataset to mean values at each site
    groups = df.groupby('site_no')
    well_info = groups.mean().copy()
    well_info = well_info.join(metadata, rsuffix='_met')
    well_info['start_dt'] = groups.datetime.min()
    well_info['end_dt'] = groups.datetime.max()

    # convert to meters; convert screen tops and botms to depths
    for col in ['well_el_m', 'head_m', 'head_std_m', 'screen_top', 'screen_botm']:
        well_info[col] *= unit_conversion
    well_info['screen_top'] = well_info['well_el_m'] - well_info['screen_top']
    well_info['screen_botm'] = well_info['well_el_m'] - well_info['screen_botm']

    # just the data, site numbers, times and aquifer
    df = df[['site_no', 'datetime', 'head_m', 'head_std_m']].copy()
    df['head_m'] *= unit_conversion
    df['head_std_m'] *= unit_conversion

    # #### trim down to only wells with both estimated water levels and standard deviation
    criteria = pd.notnull(well_info['head_m']) & pd.notnull(well_info['head_std_m'])
    well_info = well_info[criteria]

    # verify that all wells have a wellhead elevation
    assert not np.any(np.isnan(well_info.well_el_m))

    # make a histogram of open interval length
    open_interval_length = well_info.screen_top - well_info.screen_botm
    ax = open_interval_length.hist(bins=100)
    ax.set_xlabel('Well open interval length, in feet')
    ax.set_ylabel('Count')
    ax.text(.3, .7, 'median: {:.0f}\nmean: {:.0f}\nmax: {:.0f}\nmin:  {:.0f}'.format(open_interval_length.median(),
                                                                                     open_interval_length.mean(),
                                                                                     open_interval_length.max(),
                                                                                     open_interval_length.min()
                                                                                     ), transform=ax.transAxes)
    plt.savefig(out_plot)
    plt.close()

    # drop well_info with negative reported open interval
    well_info = well_info.loc[open_interval_length > 0]

    # project x, y to model crs
    x_pr, y_pr = project((well_info.lon.values, well_info.lat.values),
                          'epsg:{}'.format(src_epsg), 'epsg:{}'.format(dest_epsg))
    well_info.drop(['lon', 'lat'], axis=1, inplace=True)
    well_info['x_{}'.format(dest_epsg)], well_info['y_{}'.format(dest_epsg)] = x_pr, y_pr
    well_info['geometry'] = [Point(x, y) for x, y in zip(x_pr, y_pr)]

    # cull data to well_info in well info table
    df = df.loc[df.site_no.isin(well_info.index)].copy()

    # make unique n-character names for each observation
    # names indicate site location
    # 13 character length to allow for prefix_yyyymmm in actual observation names
    unique_obsnames = set()
    obsnames = []
    for sn in well_info.index.tolist():
        name = make_obsname(sn, unique_names=unique_obsnames,
                            maxlen=max_obsname_len)
        assert name not in unique_obsnames
        unique_obsnames.add(name)
        obsnames.append(name)
    well_info['obsprefix'] = obsnames

    # add area of interest information
    aoi_names = ['none'] * len(well_info)
    if aoi is not None:
        for name, shpname in aoi.items():
            aoi_info = shp2df(shpname)
            aoi_poly = aoi_info.geometry.values[0]
            aoi_names = [name if g.within(aoi_poly) else aoi_names[i]
                         for i, g in enumerate(well_info.geometry)]
    well_info['aoi'] = aoi_names

    # save out the results
    df2shp(well_info.drop(['x_5070', 'y_5070'], axis=1),
           out_shapefile, index=True, epsg=dest_epsg)
    print('writing {}'.format(out_info_csvfile))
    well_info.drop('geometry', axis=1).to_csv(out_info_csvfile, index=True, float_format='%.2f')
    print('writing {}'.format(out_data_csvfile))
    df.to_csv(out_data_csvfile, index=False, float_format='%.2f')   