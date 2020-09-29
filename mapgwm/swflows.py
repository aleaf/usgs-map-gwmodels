import collections
import os
from pathlib import Path
from shapely.geometry import Point, MultiPolygon
import numpy as np
import pandas as pd
from gisutils import shp2df, df2shp, project, get_values_at_points
from mfsetup.obs import make_obsname
from mfsetup.units import convert_volume_units, convert_time_units
from mapgwm.utils import makedirs


def format_site_ids(iterable, add_leading_zeros=False):
    """Cast site ids to strings"""
    str_ids = []
    for id in iterable:
        if add_leading_zeros:
            str_ids.append(format_usgs_sw_site_id(id))
        else:
            str_ids.append(str(id))
    return str_ids


def format_usgs_sw_site_id(stationID):
    """Add leading zeros to NWIS surface water sites, if they are missing.
    See https://help.waterdata.usgs.gov/faq/sites/do-station-numbers-have-any-particular-meaning.
    Zeros are only added to numeric site numbers less than 15 characters in length.
    """
    if not str(stationID).startswith('0') and str(stationID).isdigit() and \
            0 < int(str(stationID)[0]) < 10 and len(str(stationID)) < 15:
        return '0{}'.format(stationID)
    return str(stationID)


def preprocess_flows(data, metadata=None, flow_data_columns=['flow'],
                     start_date=None, active_area=None,
                     active_area_id_column=None,
                     active_area_feature_id=None,
                     source_crs=4269, dest_crs=5070,
                     datetime_col='datetime',
                     site_no_col='site_no',
                     line_id_col='line_id',
                     x_coord_col='x',
                     y_coord_col='y',
                     name_col='name',
                     flow_qualifier_column=None,
                     default_qualifier='measured',
                     include_sites=None,
                     include_line_ids=None,
                     source_volume_units='ft3',
                     source_time_units='s',
                     dest_volume_units='m3',
                     dest_time_units='d',
                     geographic_groups=None,
                     geographic_groups_col=None,
                     max_obsname_len=None,
                     add_leading_zeros_to_sw_site_nos=False,
                     column_renames=None,
                     outfile='../source_data/observations/flux_obs/preprocessed_flux_obs.csv',
                     ):
    """Preprocess stream flow observation data, for example, from NWIS or another data source that
    outputs time series in CSV format with site locations and identifiers.

    * Data are reprojected from a `source_crs` (Coordinate reference system; assumed to be in geographic coordinates)
      to the CRS of the model (`dest_crs`)
    * Data are culled to a `start_date` and optionally, a polygon or set of polygons defining the model area
    * length and time units are converted to those of the groundwater model.
    * Prefixes for observation names (with an optional length limit) that identify the location are generated
    * Preliminary observation groups can also be assigned, based on geographic areas defined by polygons
      (`geographic_groups` parameter)

    Parameters
    ----------
    data : csv file or DataFrame
        Time series of stream flow observations.
        Columns:

        ===================== ======================================
        site_no               site identifier
        datetime              measurement dates/times
        x                     x-coordinate of site
        y                     y-coordinate of site
        flow_data_columns     Columns of observed streamflow values
        flow_qualifier_column Optional column with qualifiers for flow values
        ===================== ======================================

        Notes:

        * x and y columns can alternatively be in the metadata table
        * flow_data_columns are denoted in `flow_data_columns`; multiple
          columns can be included to process base flow and total flow, or
          other statistics in tandem
        * For example, `flow_qualifier_column` may have "estimated" or "measured"
          flags denoting whether streamflows were derived from measured values
          or statistical estimates.

    metadata : csv file or DataFrame
        Stream flow observation site information.

        May include columns:

        ================= ================================================================================
        site_no           site identifier
        x                 x-coordinate of site
        y                 y-coordinate of site
        name              name of site
        line_id_col       Identifier for a line in a hydrography dataset that the site is associated with.
        ================= ================================================================================

        Notes:

        * other columns in metadata will be passed through to the metadata output

    flow_data_columns : list of strings
        Columns in data with flow values or their statistics.
        By default, ['q_cfs']
        start_date : str (YYYY-mm-dd)
        Simulation start date (cull observations before this date)
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
    datetime_col : str, optional
        Column name in data with observation date/times,
        by default 'datetime'
    site_no_col : str, optional
        Column name in data and metadata with site identifiers,
        by default 'site_no'
    line_id_col : str, optional
        Column name in data or metadata with identifiers for
        hydrography lines associated with observation sites.
        by default 'line_id'
    x_coord_col : str, optional
        Column name in data or metadata with x-coordinates,
        by default 'x'
    y_coord_col : str, optional
        Column name in data or metadata with y-coordinates,
        by default 'y'
    name_col : str, optional
        Column name in data or metadata with observation site names,
        by default 'name'
    flow_qualifier_column : str, optional
        Column name in data with flow observation qualifiers, such
        as "measured" or "estimated"
        by default 'category'
    default_qualifier : str, optional
        Default qualifier to populate flow_qualifier_column if it
        is None. By default, "measured"
    include_sites : list-like, optional
        Exclude output to these sites.
        by default, None (include all sites)
    include_line_ids : list-like, optional
        Exclude output to these sites, represented by line identifiers.
        by default, None (include all sites)
    source_volume_units : str, 'm3', 'cubic meters', 'ft3', etc.
        Volume units of the source data. By default, 'ft3'
    source_time_units : str, 's', 'seconds', 'days', etc.
        Time units of the source data. By default, 's'
    dest_volume_units : str, 'm3', 'cubic meters', 'ft3', etc.
        Volume units of the output (model). By default, 'm'
    dest_time_units : str, 's', 'seconds', 'days', etc.
        Time units of the output (model). By default, 'm'
    geographic_groups : file, dict or list-like
        Option to group observations by area(s) of interest. Can
        be a shapefile, list of shapefiles, or dictionary of shapely polygons.
        A 'group' column will be created in the metadata, and observation
        sites within each polygon will be assigned the group name
        associated with that polygon.

        For example::

            geographic_groups='../source_data/extents/CompositeHydrographArea.shp'
            geographic_groups=['../source_data/extents/CompositeHydrographArea.shp']
            geographic_groups={'cha': <shapely Polygon>}

        Where 'cha' is an observation group name for observations located within the
        the area defined by CompositeHydrographArea.shp. For shapefiles,
        group names are provided in a `geographic_groups_col`.

    geographic_groups_col : str
        Field name in the `geographic_groups` shapefile(s) containing the
        observation group names associated with each polygon.
    max_obsname_len : int or None
        Maximum length for observation name prefix. Default of 13
        allows for a PEST obsnme of 20 characters or less with
        <prefix>_yyyydd or <prefix>_<per>d<per>
        (e.g. <prefix>_2d1 for a difference between stress periods 2 and 1)
        If None, observation names will not be truncated. PEST++ does not have
        a limit on observation name length.
    add_leading_zeros_to_sw_site_nos : bool
        Whether or not to pad site numbers using the
        :func:~`mapgwm.swflows.format_usgs_sw_site_id` function.
        By default, False.
    column_renames : dict, optional
        Option to rename columns in the data or metadata that are different than those listed above.
        For example, if the data file has a 'SITE_NO' column instead of 'SITE_BADGE'::

            column_renames={'SITE_NO': 'site_no'}

        by default None, in which case the renames listed above will be used.
        Note that the renames must be the same as those listed above for
        :func:`mapgwm.swflows.preprocess_flows` to work.
    outfile : str
        Where output file will be written. Metadata are written to a file
        with the same name, with an additional "_info" suffix prior to
        the file extension.

    Returns
    -------
    data : DataFrame
        Preprocessed time series
    metadata : DataFrame
        Preprocessed metadata

    References
    ----------
    `The PEST++ Manual <https://github.com/usgs/pestpp/tree/master/documentation>`

    Notes
    -----

    """
    # outputs
    outpath, filename = os.path.split(outfile)
    makedirs(outpath)
    outname, ext = os.path.splitext(outfile)
    out_info_csvfile = outname + '_info.csv'
    out_data_csvfile = outfile
    out_shapefile = outname + '_info.shp'

    # read the source data
    if isinstance(data, str) or isinstance(data, Path):
        df = pd.read_csv(data)
    elif isinstance(data, pd.DataFrame):
        pass
    else:
        raise TypeError('Data type not understood. '
                         '"data" must be a csv file or DataFrame.')
    # check the columns
    for col in [datetime_col] + flow_data_columns:
        assert col in df.columns, "Column {} not found in {}".format(col,
                                                                     data)
    assert any({site_no_col, line_id_col}.intersection(df.columns)), \
        "Neither {} or {} found in {}. Need to specify a site_no_col or line_id_col".format(site_no_col,
                                                                                            line_id_col, data)
    # rename input columns to these names,
    # for consistent output
    dest_columns = {datetime_col: 'datetime',
                    site_no_col: 'site_no',
                    line_id_col: 'line_id',
                    x_coord_col: 'x',
                    y_coord_col: 'y',
                    name_col: 'name',
                    flow_qualifier_column: 'category'
                    }
    # update the default column renames
    # with any supplied via column_renames parameter
    if isinstance(column_renames, collections.Mapping):
        dest_columns.update(column_renames)
    df.rename(columns=dest_columns, inplace=True)
    flow_data_columns = [c if c not in dest_columns else dest_columns[c]
                         for c in flow_data_columns]
    # convert site numbers to strings;
    # add leading 0s to any USGS sites that should have them
    if 'site_no' in df.columns:
        df['site_no'] = format_site_ids(df['site_no'], add_leading_zeros_to_sw_site_nos)
    else:
        df['site_no'] = df[line_id_col]

    # read the source data
    if metadata is not None:
        if isinstance(metadata, str) or isinstance(metadata, Path):
            md = pd.read_csv(metadata)
        elif isinstance(metadata, pd.DataFrame):
            pass
        else:
            raise TypeError('Metadata type not understood. '
                             '"metadata" must be a csv file or DataFrame.')
        if site_no_col not in md.columns or 'site_no' not in df.columns:
            raise IndexError('If metadata are supplied, both data and metadata must '
                             'have a site_no column.')
        md.rename(columns=dest_columns, inplace=True)
        md['site_no'] = format_site_ids(md['site_no'], add_leading_zeros_to_sw_site_nos)
        md.index = md['site_no']
        by_site = df.groupby('site_no')
        md['start_dt'] = pd.DataFrame(by_site['datetime'].first())
    else:
        by_site = df.groupby('site_no')
        md = pd.DataFrame(by_site['datetime'].first())
        md.columns = ['start_dt']
        md['site_no'] = md.index

    md['end_dt'] = pd.DataFrame(by_site['datetime'].last())
    md['n'] = pd.DataFrame(by_site['datetime'].count())
    md.reset_index(inplace=True, drop=True)

    # assign metadata if supplied
    for col in 'x', 'y', 'line_id', 'name':
        if col in df.columns and col not in md.columns:
            by_site_no = dict(zip(df['site_no'], df[col]))
            md[col] = [by_site_no[sn] for sn in md['site_no']]
            df.drop(col, axis=1, inplace=True)

    # index the dataframe to times;
    # truncate data before start date
    df.index = pd.to_datetime(df['datetime'])
    df.index.name = 'datetime'
    df = df.loc[start_date:].copy()

    # project x, y to model crs
    x_pr, y_pr = project((md.x.values, md.y.values), source_crs, dest_crs)
    md['x'], md['y'] = x_pr, y_pr
    md['geometry'] = [Point(x, y) for x, y in zip(x_pr, y_pr)]

    # cull data to that within the model area
    if active_area is not None:
        active_area_df = shp2df(active_area, dest_crs=dest_crs)
        if active_area_id_column is not None and active_area_feature_id is not None:
            loc = active_area_df[active_area_id_column] == active_area_feature_id
            assert any(loc), "feature {} not found!".format(active_area_feature_id)
            active_area_polygon = active_area_df.loc[loc, 'geometry']
        else:
            active_area_polygon = MultiPolygon(active_area_df.geometry.tolist())
        within = np.array([g.within(active_area_polygon) for g in md.geometry])
        if not np.all(within):
            print('Culling {} wells outside of the model area defined by {}.'
                  .format(np.sum(~within), active_area))
        md = md.loc[within]
        df_within = df.site_no.isin(md['site_no'])
        df = df.loc[df_within]

    # get the hydrography IDs corresponding to each site
    # using the included lookup table
    #if 'line_id' not in df.columns:
    #    assert line_id_lookup is not None, \
    #    "need to include line_ids in a column, or line_id_lookup dictionary mapping line_ids to site numbers"
    #    df = df.loc[df['site_no'].isin(line_id_lookup)].copy()
    #    df['line_id'] = [line_id_lookup[sn] for sn in df['site_no']]
        
    if include_sites is not None:
        df = df.loc[df.site_no.isin(include_sites)]
    if include_line_ids is not None:
        df = df.loc[df.line_id.isin(include_line_ids)]
    
    # convert units
    # ensure that flow values are numeric (may be objects if taken directly from NWIS)
    unit_conversion = (convert_volume_units(source_volume_units, dest_volume_units) /
                       convert_time_units(source_time_units, dest_time_units))
    for flow_col in flow_data_columns:
        df[flow_col] = pd.to_numeric(df[flow_col], errors='coerce') * unit_conversion
    df.dropna(subset=flow_data_columns, axis=0, inplace=True)

    # reformat qualifiers for consistent output
    # (lump to dest category columns of either estimated or measured)
    # with measured including values derived from baseflow separation or actual measurements)
    # output column name for flow qualifier column:
    dest_flow_qualifier_column = 'category'
    if flow_qualifier_column is not None:
        flow_qualifiers = {'calculated': 'measured',  # 'measured',
                           'base flow separated from measured values': 'measured',  # 'measured',
                           'measured total flow': 'measured',
                           'estimated gaged': 'estimated',
                           'estimated ungaged': 'estimated'}
        df[dest_flow_qualifier_column] = df[flow_qualifier_column].replace(flow_qualifiers)
    else:
        df['category'] = default_qualifier

    # make unique n-character prefixes (site identifiers) for each observation location
    # 13 character length allows for prefix_yyyymmm in 20 character observation names
    # (BeoPEST limit)
    unique_obsnames = set()
    obsnames = []
    for sn in md['site_no'].tolist():
        if max_obsname_len is not None:
            name = make_obsname(sn, unique_names=unique_obsnames,
                                maxlen=max_obsname_len)
            assert name not in unique_obsnames
        else:
            name = sn
        unique_obsnames.add(name)
        obsnames.append(name)
    md['obsprefix'] = obsnames

    # add area of interest information
    md['group'] = 'fluxes'
    if geographic_groups is not None:
        if isinstance(geographic_groups, dict):
            pass
        else:
            geo_group_dict = {}
            if isinstance(geographic_groups, str) or isinstance(geographic_groups, Path):
                geographic_groups = [geographic_groups]
            for item in reversed(geographic_groups):
                try:
                    group_info = shp2df(str(item), dest_crs=dest_crs)
                    groups = dict(zip(group_info[geographic_groups_col],
                                      group_info['geometry']))
                    geo_group_dict.update(groups)
                except:
                    pass
        for group_name, polygon in geo_group_dict.items():
            within = [g.within(polygon) for g in md.geometry]
            md.loc[within, 'group'] = group_name

    # data columns
    data_cols = ['site_no', 'datetime'] + flow_data_columns + ['category']
    df = df[data_cols]

    # save out the results
    df2shp(md.drop(['x', 'y'], axis=1),
           out_shapefile, crs=dest_crs)
    print('writing {}'.format(out_info_csvfile))
    md.drop('geometry', axis=1).to_csv(out_info_csvfile, index=False, float_format='%g')
    print('writing {}'.format(out_data_csvfile))
    df.to_csv(out_data_csvfile, index=False, float_format='%g')
    return df, md


def aggregrate_values_to_stress_periods(data, perioddata,
                                        datetime_col='datetime',
                                        values_col='values',
                                        id_col='id',
                                        category_col='qualifier',
                                        keep_columns=None):
    """Pandas sausage-making to take flow values at arbitrary times,
    and average them to model stress periods defined in a perioddata dataframe.
    Optionally, a category column identifying flow values as 'measured' or
    'estimated' can also be read. Measured values are used for the averages
    where available; the number of measured and estimated values contributing
    to each average are tallied in the output.

    Parameters
    ----------
    data : DataFrame
        Input data
    perioddata : DataFrame
        Stress Period start/end times. Must have columns:

        ============== ================= ========================
        per            int               MODFLOW Stress Period
        start_datetime str or datetime64 Period start time
        end_datetime   str or datetime64 Period end time
        ============== ================= ========================
        
    datetime_col : str
        Column in data for Measurement dates (str or datetime64)
    id_col: int or str
        Column in data identifying the hydrography line or measurement site
        associated with flow value. Valid identifiers corresponding to the source
        hydrography (e.g. NHDPlus COMIDs) may be needed for locating the flows
        within the stream network, unless x and y coordinates are available.
    values_col : float
        Column in data with observed or estimated values.
    category_col : str; categorical
        Column in data with 'measured' or 'estimated' flags indicating how each flow
        value was derived. If None, 'measured' is used for all flows. By default, 'qualifier'.
    site_no_col : str
        Optional column in data identifying the measurement site associated with each value,
        for example, for keeping track of measurement site numbers or names that are different
        than the hydrography identifiers in line_id_col.


    Returns
    -------
    df_per : DataFrame
        Stress period averages, and metadata describing the source of the
        averages (number of estimated vs. measured values). For each site,
        also includes averages and standard deviations for measurements
        from outside the stress periods defined in perioddata. For sites
        with no measurements within a stress period, an average of all
        other measurements is used.

    Notes
    -----
    This method is similar to mfsetup.tdis.aggregate_dataframe_to_stress_period in
    what it does (resample time series to model stress periods), and for modflow-setup, 
    has is superceded by that method (called via TransientTabularData object). Keeping this
    method here though, in case we need to use it in the future.
    Key differences between two methods:
    * this method operates on the whole timeseres of model stress periods instead of a single stress period
    * this method allows for specification of both measured and estimated values; the number of estimated and measured values
    contributing to each average are included in the output
    * duplicate sites in mfsetup.tdis.aggregate_dataframe_to_stress_period (for example,
    to handle wells located in the same model cell) can be aggregated with sum, mean, first, etc., or
    by raising an error. In this method, the duplicate site with the most measurements (vs. estimates) is retained.
    * this method fills stress periods without measurements using the mean values for all time.
    """
    # make dictionary of data keyed by site
    values = data.copy()

    # optionally keep track of auxillary data
    # (site number, name, etc.)
    if keep_columns is not None:
        auxillary_info = {}
        for col in keep_columns:
            auxillary_info[col] = dict(zip(values[id_col], values[col]))

    # create category column if there is none, to conform to logic below
    if category_col not in values.columns:
        values[category_col] = 'measured'

    values[datetime_col] = pd.to_datetime(values[datetime_col])
    values.index = pd.MultiIndex.from_tuples(list(zip(values[id_col], values[datetime_col], values[category_col])))
    values = {s: values.loc[s] for s in values.index.levels[0]}

    # Compute average base flow for each site, stress period
    # fill missing values using means
    # add columns describing source of base flow values
    # vs = []  # list of processed monthly baseflow dataframes
    vpers = []  # list of processed period baseflow dataframes
    for k, v in values.items():

        # create single set of values
        # Use measured data where available; estimates where unvailable
        # Note: unstack will fail if there are two measurement sites on the same line_id
        # (with the same dates, which would be duplicates wrt the line_id)
        # look at the fraction of values that are estimated;
        # keep the site with the most measurements
        if np.any(np.any(v.index.duplicated())):
            duplicate_sites = np.array(list(set(v.site_no)))
            pct_estimated = []
            for site_no in duplicate_sites:
                site_values = v.loc[v.site_no == site_no]
                pct_estimated.append(np.sum(site_values[category_col] == 'estimated') / len(site_values))
            inds = np.argsort(pct_estimated)
            duplicate_sites = duplicate_sites[inds]
            v = v.loc[v.site_no == duplicate_sites[0]].copy()
        v = v.drop([datetime_col, category_col], axis=1).unstack(level=1)

        if 'measured' in v[values_col].columns:
            if len(v[values_col].columns) > 1:
                v['picked'] = v[values_col]['measured']
                notmeasured = np.isnan(v[values_col]['measured'])
                v.loc[notmeasured, 'picked'] = v[values_col]['estimated'].loc[notmeasured]
                v['method'] = 'measured'
                v.loc[notmeasured, 'method'] = 'estimated'
                if np.any(notmeasured):
                    assert v.loc[notmeasured].method.unique()[0] == 'estimated'
                if np.any(~notmeasured):
                    assert v.loc[~notmeasured].method.unique()[0] == 'measured'
            else:
                v['picked'] = v[values_col]['measured']
                v['method'] = 'measured'
        elif 'estimated' in v[values_col].columns:
            v['picked'] = v[values_col]['estimated']
            v['method'] = 'estimated'
        elif 'measured total flow' in v[values_col].columns:
            v['picked'] = v[values_col]['measured total flow']  # so that tallies go in "n_measured"
            v['method'] = 'measured total flow'
        else:
            pass

        v = v[['picked', 'method']].copy()
        v.columns = v.columns.droplevel(level=1)

        # Figure out which stress period each discharge value is in
        # 999 indicates values outside of the simulation timeframe
        v['per'] = -9999
        dfs = []
        for i, r in perioddata.iterrows():
            v.loc[r.start_datetime:r.end_datetime, 'per'] = r.per
            # shave end date off dates in period
            # (otherwise first month of next period will be included)
            per_end_datetime = pd.Timestamp(r.end_datetime) - pd.Timedelta(1)
            df = v.loc[r.start_datetime:per_end_datetime].copy()
            df['per'] = r.per
            dfs.append(df)
        dfs.append(v.loc[v.per == -9999])  # keep values outside of the simulation timeframe
        v = pd.concat(dfs)

        # Tally number of measured vs. estimated values for each SP
        vper = v.groupby(['per', 'method']).count().picked.unstack()
        vper.rename(columns={'estimated': 'n_estimated',
                             'measured': 'n_measured'}, inplace=True)
        for c in ['n_estimated', 'n_measured']:
            if c not in vper.columns:
                vper[c] = 0
        # compute mean value and stdev for each stress period
        vper['Q_avg'] = v.groupby('per').mean()
        vper['Q_std'] = v.groupby('per').std()
        vper.columns.name = None

        # reindex to model times dataframe
        # (pandas will fill empty periods with nans)
        # use an outer join to retain values outside of the model time period
        vper = perioddata.join(vper, how='outer')

        # fill in the start/end times for data outside of simulation period
        vper.loc[-9999, 'start_datetime'] = v.loc[v.per == -9999].index.min()
        vper.loc[-9999, 'end_datetime'] = v.loc[v.per == -9999].index.max()
        vper.loc[-9999, 'per'] = -9999
        vper.per = vper.per.astype(int)

        # nans in n_measured and n_estimated columns mean there were no values
        # for that period. Fill with zero, convert columns to ints.
        vper.fillna({'n_measured': 0, 'n_estimated': 0}, inplace=True)
        vper['n_measured'] = vper.n_measured.astype(int)
        vper['n_estimated'] = vper.n_estimated.astype(int)

        # fill in the missing data with mean value
        vper['filled'] = np.isnan(vper.Q_avg)
        vper['Q_avg'] = vper.Q_avg.fillna(vper.Q_avg.mean())

        v[id_col] = k
        vper[id_col] = k

        # map any auxillary info back to output using site numbers
        if keep_columns is not None:
            for col in keep_columns:
                vper[col] = auxillary_info[col][k]

        # vs.append(v)
        vpers.append(vper)

    df_per = pd.concat(vpers, sort=True)
    return df_per