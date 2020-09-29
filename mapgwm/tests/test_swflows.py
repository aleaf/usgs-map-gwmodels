import os
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from mapgwm.swflows import aggregrate_values_to_stress_periods, preprocess_flows


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
def test_aggregrate_values_to_stress_periods(obsdata, times, category_col, keep_columns):
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


@pytest.mark.parametrize(('datafile,metadata_file,flow_data_cols,site_no_col,x_coord_col,y_coord_col,'
                          'flow_qualifier_column,line_id_col,line_id_lookup,include_sites,source_crs'),
                         (('swflows/MAP_baseflow_and_total_flow__at_325_supplemental_points__13Feb2020.csv',
                           None, ['predicted_total_flow', 'predicted_bf'],
                           'comid', 'X', 'Y', 'category', None, None, None, 4269),
                          ('swflows/nwis_dvs.csv', 'swflows/nwis_dv_sites.csv', ['q_cfs'], 'site_no',
                           'dec_long_va', 'dec_lat_va', None, None, None, None, 4269))
                         )
def test_preprocess_flows(test_data_path, datafile, metadata_file, flow_data_cols, site_no_col,
                          x_coord_col, y_coord_col, flow_qualifier_column,
                          line_id_col, line_id_lookup, include_sites, source_crs,
                          test_output_folder):
    datafile = Path(test_data_path, datafile)
    if metadata_file is not None:
        metadata_file = Path(test_data_path, metadata_file)

    geographic_groups = [test_data_path / 'extents/CompositeHydrographArea.shp',
                         test_data_path / 'extents/MAP_generalized_regions.shp'
                         ]
    outfile = os.path.join(test_output_folder, 'preprocessed_flows.csv')
    column_renames = {'predicted_total_flow': 'qtotal_m3d',
                      'predicted_bf': 'qbase_m3d',
                      'q_cfs': 'q_m3d',
                      'station_nm': 'name'}
    data, metadata = preprocess_flows(datafile,
                                      metadata_file,
                                      flow_data_columns=flow_data_cols,
                                      start_date='2008-01-01',
                                      active_area=os.path.join(test_data_path, 'extents/ms_delta.shp'),
                                      datetime_col='datetime',
                                      site_no_col=site_no_col,
                                      line_id_col=line_id_col,
                                      x_coord_col=x_coord_col,
                                      y_coord_col=y_coord_col,
                                      flow_qualifier_column=flow_qualifier_column,
                                      line_id_lookup=None,
                                      include_sites=None,
                                      include_line_ids=None,
                                      source_crs=source_crs,
                                      source_volume_units='ft3',
                                      source_time_units='s',
                                      dest_volume_units='m3',
                                      dest_time_units='d',
                                      geographic_groups=geographic_groups,
                                      geographic_groups_col='obsgroup',
                                      column_renames=column_renames,
                                      outfile=outfile
                                     )
    assert os.path.exists(os.path.join(test_output_folder, 'preprocessed_flows.csv'))
    assert os.path.exists(os.path.join(test_output_folder, 'preprocessed_flows_info.csv'))
    expected_data_columns = ['site_no', 'datetime'] + flow_data_cols + ['category']
    expected_data_columns = [column_renames.get(c, c) for c in expected_data_columns]
    assert np.all(data.columns == expected_data_columns)
    assert data.index.name == 'datetime'
    assert not any({'site_no', 'x', 'y', 'start_dt', 'end_dt', 'n', 'name',
                    'obsprefix', 'group'}.difference(metadata.columns))

    # check that units were converted to m3/d
    flow_col = 'q_m3d'
    if 'predicted_total_flow' in flow_data_cols:
        assert np.all(data.qtotal_m3d >= data.qbase_m3d * .8)
        flow_col = 'qbase_m3d'
    assert data[flow_col].max() > 1e6
