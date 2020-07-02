import numpy as np
import pandas as pd
import pytest
from mapgwm.swflows import aggregrate_values_to_stress_periods


@pytest.fixture()
def obsdata():
    return pd.DataFrame({'datetime': ['2001-01-01', # per 0, site 2000
                                  '2002-01-01', # per 0, site 2001
                                  '2002-02-02', # per 0, site 2001
                                  '2014-10-02', # per 1, site 2002
                                  '2015-01-01', # per 1, site 2002
                                  '2015-01-02', # per 1, site 2002
                                  '2015-01-02', # per 1, site 2000
                                  '2017-01-01', #        site 2003
                                  '2018-01-01', #        site 2003
                                  '2019-01-01'], #       site 2000
                         'flow_m3d': [1, 1, 0, 4, 2, 3, 2, 1, 0, 10],
                         'comment': ['measured',
                                     'measured',
                                     'measured',
                                     'estimated',
                                     'estimated',
                                     'estimated',
                                     'estimated',
                                     'measured',
                                     'measured',
                                     'estimated'
                                     ],
                         'site_no': [2000,
                                     2001,
                                     2001,
                                     2002,
                                     2002,
                                     2002,
                                     2000,
                                     2003,
                                     2003,
                                     2000
                                     ],
                         'line_id': [1002000,
                                     1002001,
                                     1002001,
                                     1002002,
                                     1002002,
                                     1002002,
                                     1002000,
                                     1002003,
                                     1002003,
                                     1002000
                                     ]
                         })


@pytest.fixture()
def times():
    start_datetime = pd.to_datetime(['2001-01-01', '2014-10-01'])
    end_datetime = pd.to_datetime(['2014-09-30', '2015-09-30'])
    perlen = [(edt - sdt).days for edt, sdt in zip(start_datetime, end_datetime)]
    times = pd.DataFrame({'start_datetime': start_datetime,
                          'end_datetime': end_datetime,
                          'per': [0, 1],
                          'perlen': perlen,
                          'steady': False})
    return times


@pytest.mark.parametrize('category_col', (None, 'comment'))
@pytest.mark.parametrize('keep_columns', (None, ['site_no']))
def test_process_flows(obsdata, times, category_col, keep_columns):
    results = aggregrate_values_to_stress_periods(obsdata, times,
                                                  datetime_col='datetime',
                                                  values_col='flow_m3d',
                                                  id_col='line_id',
                                                  category_col=category_col,
                                                  keep_columns=keep_columns
                                                  )
    sortby = ['per']
    if keep_columns is not None and 'site_no' in keep_columns:
        sortby.append('site_no')
    results.sort_values(by=sortby, inplace=True)
    # there should be results for each site, for each period
    # including -9999 period that represents anything outside of the
    # model timeframe
    if keep_columns is not None and 'site_no' in keep_columns:
        assert np.array_equal(results.site_no.values,
                              np.array([2000, 2001, 2002, 2003]*3))
    # rows without measurements should be indicated as filled
    assert np.array_equal(((results.n_estimated == 0) & (results.n_measured == 0)),
                          results.filled.values)
    # expected values for filled based on whether site had any measurements in period
    assert results.filled.tolist() == [False, True, True, False, # 2000 and 2003 both had measurements outside model period
                                       False, False, True, True,
                                       False, True, False, True]
    assert np.array_equal(results.Q_avg.values,
                          np.array([10.00000, # 2000 value outside model period
                                    0.500000,
                                    3.000000,
                                    0.500000,
                                    1.000000, # 2000 value in per 0
                                    0.500000,
                                    3.000000,
                                    0.500000,
                                    2.000000, # 2000 value in per 1
                                    0.500000,
                                    3.000000,
                                    0.500000]))
    # standard deviations should be nan in rows with fewer than two measurements
    assert np.array_equal(((results.n_estimated < 2) & (results.n_measured < 2)),
                          np.isnan(results.Q_std.values))
    if keep_columns is not None and 'site_no' in keep_columns:
        assert results.loc[(results.per == -9999) &
                           (results.site_no == 2003), 'Q_std'].values[0] == \
               np.std([0, 1], ddof=1)
        assert results.loc[(results.per == 1) &
                           (results.site_no == 2002), 'Q_std'].values[0] == \
               np.std([4, 3, 2], ddof=1)