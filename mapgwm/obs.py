"""General functions for working with observations.
"""
from collections.abc import Mapping
import os
from shapely.geometry import Point
import pandas as pd
from gisutils import df2shp, project
from mfsetup.obs import make_obsname
from mfsetup.units import convert_length_units, convert_volume_units, convert_time_units
from mapgwm.utils import makedirs, assign_geographic_obsgroups, cull_data_to_active_area


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


def preprocess_obs(data, metadata=None, data_columns=['flow'],
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
                     qualifier_column=None,
                     default_qualifier='measured',
                     obstype='flow',
                     include_sites=None,
                     include_line_ids=None,
                     source_length_units='ft',
                     source_time_units='s',
                     dest_length_units='m',
                     dest_time_units='d',
                     geographic_groups=None,
                     geographic_groups_col=None,
                     max_obsname_len=None,
                     add_leading_zeros_to_sw_site_nos=False,
                     column_renames=None,
                     outfile=None,
                     ):
    """Preprocess observation data, for example, from NWIS or another data source that
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
        Time series of observations.
        Columns:

        ===================== ======================================
        site_no               site identifier
        datetime              measurement dates/times
        x                     x-coordinate of site
        y                     y-coordinate of site
        data_columns          Columns of observed values
        qualifier_column      Optional column with qualifiers for values
        ===================== ======================================

        Notes:

        * x and y columns can alternatively be in the metadata table
        * data_columns are denoted in `data_columns`; multiple
          columns can be included to process base flow and total flow, or
          other statistics in tandem
        * For example, `qualifier_column` may have "estimated" or "measured"
          flags denoting whether streamflows were derived from measured values
          or statistical estimates.

    metadata : csv file or DataFrame
        Observation site information.

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

    data_columns : list of strings
        Columns in data with values or their statistics.
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
    qualifier_column : str, optional
        Column name in data with observation qualifiers, such
        as "measured" or "estimated"
        by default 'category'
    default_qualifier : str, optional
        Default qualifier to populate qualifier_column if it
        is None. By default, "measured"
    obstype : str, optional
        Modflow-6 observation type (e.g. 'downstream-flow' or 'stage'). 
        The last part of the name (after the last hyphen) is used as a suffix in the output 
        ``obsprefix`` column. E.g. 07275000-flow for downstream or upstream-flow at site 07275000.
        By default, 'flow'
    include_sites : list-like, optional
        Exclude output to these sites.
        by default, None (include all sites)
    include_line_ids : list-like, optional
        Exclude output to these sites, represented by line identifiers.
        by default, None (include all sites)
    source_length_units : str, 'm3', 'm', 'cubic meters', 'ft3', etc.
        Length or volume units of the source data. By default, 'ft3'
    source_time_units : str, 's', 'seconds', 'days', etc.
        Time units of the source data. By default, 's'
    dest_length_units : str, 'm3', 'cubic meters', 'ft3', etc.
        Length or volume units of the output (model). By default, 'm'
    dest_time_units : str, 's', 'seconds', 'days', etc.
        Time units of the output (model). By default, 'd'
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
        :func:`mapgwm.swflows.preprocess_obs` to work.
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
    if outfile is not None:
        outpath, filename = os.path.split(outfile)
        makedirs(outpath)
        outname, ext = os.path.splitext(outfile)
        out_info_csvfile = outname + '_info.csv'
        out_data_csvfile = outfile
        out_shapefile = outname + '_info.shp'

    # read the source data
    if not isinstance(data, pd.DataFrame):
        df = pd.read_csv(data, dtype={site_no_col: object})
    else:
        df = data.copy()
    # check the columns
    for col in [datetime_col] + data_columns:
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
                    qualifier_column: 'category'
                    }
    # update the default column renames
    # with any supplied via column_renames parameter
    if isinstance(column_renames, Mapping):
        dest_columns.update(column_renames)
    df.rename(columns=dest_columns, inplace=True)
    data_columns = [c if c not in dest_columns else dest_columns[c]
                         for c in data_columns]
    # convert site numbers to strings;
    # add leading 0s to any USGS sites that should have them
    if 'site_no' in df.columns:
        df['site_no'] = format_site_ids(df['site_no'], add_leading_zeros_to_sw_site_nos)
    else:
        df['site_no'] = df[line_id_col]
    
    # make obsprefix names with site and observation type
    df['obsprefix'] = [f"{site_no}-{obstype.split('-')[-1]}" 
                       for site_no in df['site_no']]

    # read the source data
    if metadata is not None:
        if not isinstance(metadata, pd.DataFrame):
            md = pd.read_csv(metadata, dtype={site_no_col: object})
        else:
            md = metadata.copy()
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
            if col != 'line_id':
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
        df, md = cull_data_to_active_area(df, active_area,
                                          active_area_id_column,
                                          active_area_feature_id,
                                          data_crs=dest_crs, metadata=md)

    # get the hydrography IDs corresponding to each site
    # using the included lookup table
    #if 'line_id' not in df.columns:
    #    assert line_id_lookup is not None, \
    #    "need to include line_ids in a column, or line_id_lookup dictionary mapping line_ids to site numbers"
    #    df = df.loc[df['site_no'].isin(line_id_lookup)].copy()
    #    df['line_id'] = [line_id_lookup[sn] for sn in df['site_no']]
        
    if include_sites is not None:
        md_included = md.site_no.isin(include_sites)
        if not any(md_included):
            raise ValueError('None of the sites in include_sites are in the metadata or model area!')
        md = md.loc[md_included]
        data_included = df.site_no.isin(include_sites)
        if not any(data_included):
            raise ValueError('None of the sites in include_sites are in the data or model area!')
        df = df.loc[data_included]
    if include_line_ids is not None:
        md_included = md.line_id.isin(include_line_ids)
        if not any(md_included):
            raise ValueError('None of the sites in include_sites are in the metadata or model area!')
        md = md.loc[md_included]
        data_included = df.line_id.isin(include_line_ids)
        if not any(data_included):
            raise ValueError('None of the sites in include_sites are in the data or model area!')
        df = df.loc[data_included]
    
    # convert units
    # ensure that values are numeric (may be objects if taken directly from NWIS)
    if obstype == 'stage':
        unit_conversion = convert_length_units(source_length_units, dest_length_units)
    else:
        unit_conversion = (convert_volume_units(source_length_units, dest_length_units) /
                        convert_time_units(source_time_units, dest_time_units))
    for obs_col in data_columns:
        df[obs_col] = pd.to_numeric(df[obs_col], errors='coerce') * unit_conversion
    df.dropna(subset=data_columns, axis=0, inplace=True)

    # reformat qualifiers for consistent output
    # (lump to dest category columns of either estimated or measured)
    # with measured including values derived from baseflow separation or actual measurements)
    # output column name for qualifier column:
    dest_qualifier_column = 'category'
    if qualifier_column is not None:
        qualifiers = {'calculated': 'measured',  # 'measured',
                           'base flow separated from measured values': 'measured',  # 'measured',
                           'measured total flow': 'measured',
                           'estimated gaged': 'estimated',
                           'estimated ungaged': 'estimated'}
        df[dest_qualifier_column] = df[qualifier_column].replace(qualifiers)
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
        name = name + f"-{obstype.split('-')[-1]}"
        unique_obsnames.add(name)
        obsnames.append(name)
    md['obsprefix'] = obsnames

    # add area of interest information
    md['group'] = 'fluxes'
    md = assign_geographic_obsgroups(md, geographic_groups,
                                     geographic_groups_col,
                                     metadata_crs=dest_crs)

    # data columns
    data_cols = ['site_no', 'line_id', 'datetime', 'obsprefix'] + data_columns + ['category']
    #if 'line_id' in md.columns and 'line_id' not in df.columns:
    #    # only map line_ids to data if there are more site numbers
    #    # implying that no site number maps to more than one line_id
    #    if len(set(df.site_no)) >= len(set(df.line_id)):
    #        ids = dict(zip(md['site_no'], md['line_id']))
    #    df['line_id'] = [ids[sn] for sn in df['site_no']]
    data_cols = [c for c in data_cols if c in df.columns]
    df = df[data_cols]

    md.index = md['site_no']
    # save out the results
    if outfile is not None:
        df2shp(md.drop(['x', 'y'], axis=1),
               out_shapefile, crs=dest_crs)
        print('writing {}'.format(out_info_csvfile))
        md.drop('geometry', axis=1).to_csv(out_info_csvfile, index=False, float_format='%g')
        print('writing {}'.format(out_data_csvfile))
        df.to_csv(out_data_csvfile, index=False, float_format='%g')
    return df, md

