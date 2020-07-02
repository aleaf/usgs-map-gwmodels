import numpy as np
import pandas as pd
from mfsetup.units import convert_volume_units, convert_time_units


def str_ids(iterable):
    """Cast site ids to strings"""
    str_ids = []
    for id in iterable:
        if not str(id).startswith('0') and str(id).isdigit():
            str_ids.append(correct_station_id(id))
        else:
            str_ids.append(str(id))
    return str_ids


def correct_station_id(stationID):
    """For dealing with USGS site IDs that have leading 0s.
    """
    if 1 < int(str(stationID)[0]) < 10 and len(str(stationID)) < 15:
        return '0{}'.format(stationID)
    return str(stationID)


def preprocess_flows(streamflows_file,
                     start_date=None,
                     datetime_col='datetime',
                     site_no_col='site_no',
                     line_id_col='line_id',
                     flows_col='flow',
                     flow_qualifier_column=None,
                     line_id_lookup=None,
                     include_sites=None,
                     include_line_ids=None,
                     source_volume_units='ft3',
                     source_time_units='s',
                     dest_volume_units='m3',
                     dest_time_units='d'
                     ):
    
    # read the source data
    df = pd.read_csv(streamflows_file)
    # check the columns
    for col in datetime_col, flows_col:
        assert col in df.columns, "Column {} not found in {}".format(col, 
                                                                     streamflows_file)
    assert any({site_no_col, line_id_col}.intersection(df.columns)), \
        "Neither {} or {} found in {}. Need to specify a site_no_col or line_id_col".format(site_no_col, line_id_col, streamflows_file)
    
    # rename input columns to these names,
    # for consistent output
    dest_columns = {datetime_col: 'datetime',
                    site_no_col: 'site_no',
                    line_id_col: 'line_id',
                    flows_col: 'flow_{}'.format(source_volume_units)
                    }
    df.rename(columns=dest_columns, inplace=True)
    
    # convert site numbers to strings;
    # add leading 0s to any USGS sites that should have them
    if 'site_no' in df.columns:
        df['site_no'] = str_ids(df['site_no'])
    
    # index the dataframe to times;
    # truncate data before start date
    df.index = pd.to_datetime(df['datetime'])
    df.index.name = 'datetime'
    df = df.loc[start_date:].copy()

    # get the hydrography IDs corresponding to each site
    # using the included lookup table
    if 'line_id' not in df.columns:
        assert line_id_lookup is not None, \
        "need to include line_ids in a column, or line_id_lookup dictionary mapping line_ids to site numbers"
        df = df.loc[df['site_no'].isin(line_id_lookup)].copy()
        df['line_id'] = [line_id_lookup[sn] for sn in df['site_no']]
        
    if include_sites is not None:
        df = df.loc[df.site_no.isin(include_sites)]
    if include_line_ids is not None:
        df = df.loc[df.line_id.isin(include_line_ids)]
    
    # convert units
    unit_conversion = (convert_volume_units(source_volume_units, dest_volume_units) /
                       convert_time_units(source_time_units, dest_time_units))
    converted_flows_col = 'flow_' + dest_volume_units + dest_time_units
    df[converted_flows_col] = df['flow_{}'.format(source_volume_units)] * unit_conversion

    # reformat qualifiers for consistent output
    # (lump to either estimated or measured, 
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
    
    output_columns = [c for c in ('site_no', 'line_id', converted_flows_col, dest_flow_qualifier_column)
                      if c in df.columns]
    return df[output_columns]


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
    perioddata : DataFrame
        Stress Period start/end times. Must have columns:
        per : int
            MODFLOW Stress Period
        start_datetime : str or datetime64
            Period start time
        end_datetime : str or datetime64
            Period end time
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