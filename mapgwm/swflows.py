"""
Code for preprocessing streamflow data.
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from mapgwm.obs import preprocess_obs


def preprocess_flows(*args, **kwargs):
    warnings.warn('The preprocess_flows() function is deprecated; use mapgwm.obs.preprocess_obs() instead.')
    return preprocess_obs(*args, **kwargs)


def combine_measured_estimated_values(measured_values, estimated_values,
                                      measured_values_data_col, estimated_values_data_col,
                                      dest_values_col='obsval',
                                      resample_freq='MS', how='mean'):
    """Combine time series of measured and estimated values for multiple sites,
    giving preference to measured values.

    Parameters
    ----------
    measured_values : csv file or DataFrame
        Time series of measured values at multiple sites, similar to that
        output by :func:`~mapgwm.swflows.preprocess_flows`.

        Columns:

        ===================== ===========================================
        site_no               site identifiers; read-in as strings
        datetime              measurement dates/times
        data columns          columns with floating-point data to combine
        ===================== ===========================================

    estimated_values : csv file or DataFrame
        Time series of measured values at multiple sites, similar to that
        output by :func:`~mapgwm.swflows.preprocess_flows`.

        Columns:

        ===================== ===========================================
        site_no               site identifiers; read-in as strings
        datetime              measurement dates/times
        data columns          columns with floating-point data to combine
        ===================== ===========================================

    measured_values_data_col : str
        Column in `measured_values` with data to combine.
    estimated_values_data_col : str
        Column in `estimated_values` with data to combine.
    dest_values_col : str
        Output column with combined data from `measured_values_data_col`
        and estimated_values_data_col, by default 'obsval'
    resample_freq : str or DateOffset
        Any `pandas frequency alias <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timeseries-offset-aliases>`_
        The data columns in `measured_values` and `estimated_values` are resampled
        to this fequency using the method specified by ``how``.
        By default, 'MS' (month-start)
    how : str
        Resample method. Can be any of the method calls on the pandas
        `Resampler <https://pandas.pydata.org/pandas-docs/stable/reference/resampling.html>`_
        object. By default, 'mean'

    Returns
    -------
    combined : DataFrame
        DataFrame containing all columns from `estimated_values`, the data columns
        from `measured_values`, and a `dest_values_col` consisting of measured values where
        present, and estimated values otherwise. An ``"est_"`` prefix is added to the
        estimated data columns, and a ``"meas"`` prefix is added to the measured data columns.

        Example:

        ======== ========== ========= ============== ============== ============
        site_no  datetime   category  est_qbase_m3d  meas_qbase_m3d obsval
        ======== ========== ========= ============== ============== ============
        07288000 2017-10-01 measured  47872.1        28438.7        28438.7
        07288000 2017-11-01 measured  47675.9        24484.5        24484.5
        ======== ========== ========= ============== ============== ============

        Where ``category`` denotes whether the value in obsval is measured or estimated.

    Notes
    -----
    All columns with a floating-point dtype are identified as "Data columns," and
    are resampled as specified by the `resample_freq` and `how` arguments. For all other
    columns, the first value for each time at each site is used. The resampled measured
    data are joined to the resampled estimated data on the basis of site numbers
    and times (as a pandas `MultiIndex`).


    """
    # read in the data
    if not isinstance(measured_values, pd.DataFrame):
        measured = pd.read_csv(measured_values, dtype={'site_no': object})
        measured['datetime'] = pd.to_datetime(measured['datetime'])
        measured.index = measured['datetime']
    else:
        measured = measured_values.copy()
    if not isinstance(estimated_values, pd.DataFrame):
        df = pd.read_csv(estimated_values, dtype={'site_no': object})
        df['datetime'] = pd.to_datetime(df['datetime'])
        df.index = df['datetime']
    else:
        df = estimated_values.copy()

    # resample both timeseries to the resample_freq
    df_rs = resample_group_timeseries(df, resample_freq=resample_freq, how=how,
                                      add_data_prefix='est')
    df_rs['category'] = 'estimated'
    measured_rs = resample_group_timeseries(measured, resample_freq=resample_freq, how=how,
                                            add_data_prefix='meas')
    # fill any nan values created in resampling with 'measured' classifier
    measured_rs['category'] = 'measured'

    # add the measured values to the estimated
    measured_data_columns = measured_rs.select_dtypes(include=[np.float]).columns
    # join the data columns
    combined = df_rs.join(measured_rs[measured_data_columns], how='outer')
    # update the site_no, datetime and category columns from measured
    # (that were not filled in outer join)
    combined.update(measured_rs)

    # make the obsval column, starting with the estimated values
    combined[dest_values_col] = combined['est_' + estimated_values_data_col]
    # prefix was added to measured data column by resample_group_timeseries
    measured_values_data_col = 'meas_' + measured_values_data_col
    # populate obsval column with all measured values
    # (over-writing any pre-existing estimated values)
    has_measured = ~combined[measured_values_data_col].isna()
    combined.loc[has_measured, dest_values_col] = combined.loc[has_measured, measured_values_data_col]

    # verify that only nan values in obsval column are for measured sites
    # (that don't have estimated values)
    assert np.all(combined.loc[combined.obsval.isna(), 'category'] == 'measured')
    # drop observations with no measured or estimated values
    combined.dropna(subset=['obsval'], axis=0, inplace=True)
    return combined.reset_index(drop=True)


def resample_group_timeseries(df, resample_freq='MS', how='mean',
                              add_data_prefix=None):
    """Resample a DataFrame with both groups (e.g. measurement sites)
    and time series (measurements at each site).

    Parameters
    ----------
    df : DataFrame
        Time series of values at multiple sites, similar to that
        output by :func:`~mapgwm.swflows.preprocess_flows`.

        Columns:

        ===================== ===========================================
        site_no               site identifiers; read-in as strings
        datetime              measurement dates/times
        data columns          columns with floating-point data to resample
        ===================== ===========================================

    resample_freq : str or DateOffset
        Any `pandas frequency alias <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#timeseries-offset-aliases>`_
        The data columns in `measured_values` and `estimated_values` are resampled
        to this fequency using the method specified by ``how``.
        By default, 'MS' (month-start)
    how : str
        Resample method. Can be any of the method calls on the pandas
        `Resampler <https://pandas.pydata.org/pandas-docs/stable/reference/resampling.html>`_
        object. By default, 'mean'
    add_data_prefix : str
        Option to add prefix to data columns. By default, None

    Returns
    -------
    resampled : DataFrame
        Resampled data at each site.

    Notes
    -----
    All columns with a floating-point dtype are identified as "Data columns," and
    are resampled as specified by the `resample_freq` and `how` arguments. For all other
    columns, the first value for each time at each site is used.
    """
    if not isinstance(df.index, pd.DatetimeIndex):
        raise ValueError("Input must have a DatetimeIndex.")

    df_resampler = df.groupby('site_no').resample(resample_freq)
    df_rs = df_resampler.first().copy()
    # resample the data columns
    data_columns = df_rs.select_dtypes(include=[np.float]).columns
    non_data_columns = [c for c in df.columns if c not in data_columns]
    resampled_values = getattr(df_resampler[data_columns], how)()
    # put the non-data and data columns back together
    df_rs = df_rs[non_data_columns].join(resampled_values)
    # add 'est' suffix to column names
    if add_data_prefix is not None:
        for col in data_columns:
            df_rs[add_data_prefix + '_' + col] = df_rs[col]
            df_rs.drop(col, axis=1, inplace=True)
    # index dates may vary depending on freq (e.g. 'MS' vs 'M')
    df_rs['site_no'] = df_rs.index.get_level_values(0)
    df_rs['datetime'] = df_rs.index.get_level_values(1)
    return df_rs


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