from pathlib import Path

import pandas as pd
from shapely.geometry import Point

from gisutils import df2shp, project, get_values_at_points
from mfsetup.units import convert_length_units, convert_volume_units
from mapgwm.utils import cull_data_to_active_area


def read_te_water_use_spreadsheet(xlsx_file, date='2010', site_no_col='site_no',
                                  site_name_col='site_name',
                                  x_coord_col='x', y_coord_col='y',
                                  q_col='q', source_name_col='source_name',
                                  source_code_col='WATER_SOURCE_CODE',
                                  source_code_filter=['GW'],
                                  **kwargs):
    """Read water use data for thermoelectric power generation from
    a spreadsheet; filter by a source code, and rename columns to
    uniform names so that multiple datasets from spreadsheets with different
    formatting can be easily combined.

    Parameters
    ----------
    xlsx_file : Excel spreadsheet
        Thermoelectric power water use data, in a format readable by
        :func:`pandas.read_excel`.
    date : str
        Date string in a format readable by :func:`pandas.Timestamp`,
        indicating the time of the water use data (most likely a year).
    site_no_col : str
        Column in `xlsx_file` with identifiers for water users (power plants).
    site_name_col : str
        Column in `xlsx_file` with names of water users (power plants).
    x_coord_col : str
        Column in `xlsx_file` with power plant location x-coordinates.
    y_coord_col : str
        Column in `xlsx_file` with power plant location y-coordinates.
    q_col : str
        Column in `xlsx_file` with power plant pumping rates.
    source_name_col : str
        Column in `xlsx_file` with names of water sources.
    source_code_col : str
        Column in `xlsx_file` with codes for water supply sources.
    source_code_filter : sequence of strings
        Source code(s) in `source_code_col` to filter on; e.g. 'GW' for groundwater
    kwargs : keywword arguments to :func:`pandas.read_excel`
        Arguments to read_excel that might be needed to read
        `xlsx_file`; for example ``sheet_name`` or ``skiprows``.

    Returns
    -------
    df : DataFrame
        TE data with unified column names:

        ========= ====================================================
        site_no   power plant identifier (plant code)
        date      pandas datetime representative of flux (e.g. '2010')
        x         x-coordinate of withdrawl, in `source_crs`
        y         y-coordinate of withdrawl, in `source_crs`
        q         withdrawl flux, in `data_volume_units` per days
        site_name name of power plant, if provided
        ========= ====================================================

    """
    df = pd.read_excel(xlsx_file, **kwargs)
    renames = {site_no_col: 'site_no',
               site_name_col: 'site_name',
               source_name_col: 'source_name',
               source_code_col: 'source_cd',
               x_coord_col: 'x',
               y_coord_col: 'y',
               q_col: 'q'}
    df.rename(columns=renames, inplace=True)

    # filter by source code
    criteria = (df['source_cd'].isin(source_code_filter))
    if source_code_filter == 'GW':
        criteria = criteria | (df['source_name'].str.contains('well'))
    df = df.loc[criteria]

    df['date'] = pd.to_datetime(date)
    columns = ['site_no', 'date', 'x', 'y', 'q', 'site_name', 'source_name']
    columns = [c for c in columns if c in df.columns]
    return df[columns]


def preprocess_te_wateruse(data,
                           start_date=None,
                           end_date=None,
                           active_area=None,
                           active_area_id_column=None,
                           active_area_feature_id=None,
                           estimated_production_zone_top=None,
                           estimated_production_zone_botm=None,
                           estimated_production_surface_units='feet',
                           source_crs=4269,
                           dest_crs=5070,
                           interp_method='linear',
                           data_volume_units='mgal',
                           model_length_units='meters',
                           outfile=None):
    """Preprocess water use data from thermoelectric power plants:

    * reproject data to a destination CRS `dest_crs`)
    * cull data to an area of interest (`active_area`)
    * if input data do not have information on the well screen intervals;
      sample screen tops and bottoms from raster surfaces bounding
      an estimated production zone (e.g. `estimated_production_zone_top`)
    * reindex the data to continous monthly values extending from `start_date`
      to `end_date`. Typically, these would bracket the time period for which
      the pumping should be simulated in a model. For example, the earliest data
      may be from 2010, but if the model starts in 2008, it may be appropriate to
      begin using the 2010 rates then (``start_date='2008'``). If no start or end
      date are given, the first and last years of pumping in `data` are used.
    * fill empty months by interpolation via a specified `interp_method`
    * backfill any remaining empty months going back to the `start_date`
    * write processed data to a CSV file and shapefile of the same name

    Parameters
    ----------
    data : DataFrame
        Thermoelectric water use data in the following format
        (similar to that output by :func:`mapgwm.te_wateruse.read_te_water_use_spreadsheet`):

        ======= ====================================================
        site_no power plant identifier (plant code)
        date    pandas datetime representative of flux (e.g. '2010')
        x       x-coordinate of withdrawl, in `source_crs`
        y       y-coordinate of withdrawl, in `source_crs`
        q       withdrawl flux, in `data_volume_units` per days
        ======= ====================================================

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
    estimated_production_zone_top : file path
        Raster surface for assigning screen tops
    estimated_production_zone_botm : file path
        Raster surface for assigning screen bottoms
    estimated_production_surface_units : str, {'meters', 'ft', etc.}
        Length units of elevations in estimated production surface rasters.
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
    interp_method : str
        Interpolation method to use for filling pumping rates to monthly values.
        By default, 'linear'
    data_volume_units : str; e.g. 'mgal', 'm3', 'cubic feet', etc.
        Volume units of pumping data. All time units are assumed to be in days.
    model_length_units : str; e.g. 'feet', 'm', 'meters', etc.
        Length units of model.
    outfile : str
        Path for output file. A shapefile of the same name is also written.
        If None, no output file is written. By default, None

    Returns
    -------
    df_monthly : DataFrame
        

    Notes
    -----
    * time units for TE data and model are assumed to be days

    """
    df = data.copy()

    # reproject to dest_crs
    x, y = project(zip(df['x'], df['y']), source_crs, dest_crs)
    df['x'], df['y'] = x, y
    df['geometry'] = [Point(x, y) for x, y in zip(x, y)]

    # drop wells with no location information (for now)
    df.dropna(subset=['x', 'y'], axis=0, inplace=True)

    # cull sites to those within the Delta footprint
    # cull data to that within the model area
    if active_area is not None:
        df = cull_data_to_active_area(df, active_area,
                                      active_area_id_column,
                                      active_area_feature_id,
                                      data_crs=dest_crs)

    # get top and bottom of estimated production interval at each well
    if estimated_production_zone_top is not None and \
            estimated_production_zone_botm is not None:
        surf_unit_conversion = convert_length_units(estimated_production_surface_units,
                                                    model_length_units)
        x, y = df.x.values, df.y.values
        est_screen_top = get_values_at_points(estimated_production_zone_top, x, y,
                                              points_crs=dest_crs)
        est_screen_top *= surf_unit_conversion
        est_screen_botm = get_values_at_points(estimated_production_zone_botm, x, y,
                                               points_crs=dest_crs)
        est_screen_botm *= surf_unit_conversion
        df['screen_top'] = est_screen_top
        df['screen_botm'] = est_screen_botm

    # distribute fluxes to monthly values
    # set start and end dates if not already set
    if start_date is None:
        start_date = df.date.min()
    if end_date is None:
        end_date = df.date.mmax()
    groups = df.groupby('site_no')
    all_groups = []
    for site_no, group in groups:
        dfg = group.copy()

        # create a continuous monthly time index
        # labeled at the month start
        all_dates = pd.date_range(start_date, end_date, freq='MS')
        dfg.index = dfg['date']
        dfg = dfg.reindex(all_dates)

        # interpolate the discharge values;
        # back filling to the start date
        dfg['q'] = dfg.q.interpolate(method=interp_method).bfill()
        dfg['q'] *= convert_volume_units(data_volume_units, model_length_units)

        # fill remaining columns
        dfg['date'] = dfg.index
        fill_columns = set(dfg.columns).difference({'q', 'date'})
        fill_values = group.iloc[0].to_dict()
        for c in fill_columns:
            dfg[c] = fill_values[c]
        all_groups.append(dfg)
    df_monthly = pd.concat(all_groups)

    # clean up the columns
    cols = ['site_no', 'date', 'x', 'y', 'screen_top', 'screen_botm', 'q', 'geometry']
    cols += list(set(df_monthly.columns).difference(cols))
    df_monthly = df_monthly[cols]

    # write the output
    if outfile is not None:
        outfile = Path(outfile)
        df_monthly.drop('geometry', axis=1).to_csv(outfile, index=False, float_format='%g')
        print('wrote {}'.format(outfile))

        # write only unique pumping values to shapefile
        to_shapefile = df_monthly.groupby(['site_no', 'q']).first().reset_index()
        shapefile = outfile.with_suffix('.shp')
        df2shp(to_shapefile, shapefile, crs=dest_crs)
    return df_monthly
