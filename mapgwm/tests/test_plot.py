from pathlib import Path
import glob
from shapely.geometry import Point
import pytest
from mapgwm.plot import plot_wateruse
from mapgwm.utils import cull_data_to_active_area
from .test_iwum import preprocessed_iwum_data, ncfile
from .test_swuds import preprocessed_swuds_data
from .test_te_wateruse import preprocessed_te_data, te_data


def test_plot_wateruse(preprocessed_iwum_data, preprocessed_swuds_data,
                       test_data_path):

    extent = test_data_path / 'extents/shellmound_bbox.shp'
    perioddata = test_data_path / 'shellmound/tables/stress_period_data.csv'
    preprocessed_iwum_data['geometry'] = [Point(x, y) for x, y in zip(preprocessed_iwum_data.x,
                                                                      preprocessed_iwum_data.y)]
    iwum_data = cull_data_to_active_area(preprocessed_iwum_data,
                                         active_area=extent, data_crs=5070)
    swuds_data = preprocessed_swuds_data.df
    swuds_data['geometry'] = [Point(x, y) for x, y in zip(swuds_data.x, swuds_data.y)]
    swuds_data = cull_data_to_active_area(swuds_data,
                                         active_area=extent, data_crs=5070)
    wel_files = sorted(glob.glob(str(test_data_path / 'shellmound/external/wel_*.dat')))
    wel_files = dict(zip(range(1, len(wel_files) + 1), wel_files))

    add_data = {'IWUM estimates': {'data': iwum_data},
                'SWUDs data': {'data': swuds_data}
                }
    results = plot_wateruse(wel_files, perioddata, add_data)
    j=2