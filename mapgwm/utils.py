"""
Utility functions shared by other modules
"""
import os
from pathlib import Path

import gisutils
import numpy as np
from scipy import interpolate
from shapely.geometry import MultiPolygon
import pandas as pd
from gisutils import shp2df


def makedirs(path):
    dirs = [os.path.join(path, 'shps/'),
            os.path.join(path, 'figures/')]
    for folder in dirs:
        if not os.path.isdir(folder):
            os.makedirs(folder)


def assign_geographic_obsgroups(metadata, geographic_groups, geographic_groups_col,
                                metadata_crs):

    md = metadata.copy()
    if geographic_groups is not None:
        if isinstance(geographic_groups, dict):
            pass
        else:
            geo_group_dict = {}
            if isinstance(geographic_groups, str) or isinstance(geographic_groups, Path):
                geographic_groups = [geographic_groups]
            for item in reversed(geographic_groups):
                try:
                    group_info = shp2df(str(item), dest_crs=metadata_crs)
                    groups = dict(zip(group_info[geographic_groups_col],
                                      group_info['geometry']))
                    geo_group_dict.update(groups)
                except:
                    pass
        for group_name, polygon in geo_group_dict.items():
            within = [g.within(polygon) for g in md.geometry]
            md.loc[within, 'group'] = group_name
    return metadata


def cull_data_to_active_area(data, active_area, active_area_id_column,
                             active_area_feature_id,
                             data_crs, metadata=None):
    df = data.copy()
    if metadata is not None:
        md = metadata.copy()
    if isinstance(active_area, Path) or isinstance(active_area, str):
        active_area = [active_area]
    active_area = [str(filepath) for filepath in active_area]
    active_area_df = shp2df(active_area, dest_crs=data_crs)
    if active_area_id_column is not None and active_area_feature_id is not None:
        loc = active_area_df[active_area_id_column] == active_area_feature_id
        assert any(loc), "feature {} not found!".format(active_area_feature_id)
        active_area_polygon = active_area_df.loc[loc, 'geometry']
    else:
        active_area_polygon = MultiPolygon(active_area_df.geometry.tolist())

    if metadata is not None:
        within = np.array([g.within(active_area_polygon) for g in md.geometry])
        md = md.loc[within]
        df_within = df.site_no.isin(md['site_no'])
    else:
        within = np.array([g.within(active_area_polygon) for g in df.geometry])
        df_within = within
    if not np.all(within):
        print('Culling {} wells outside of the model area defined by {}.'
              .format(np.sum(~within), active_area))
    df = df.loc[df_within]
    if metadata is not None:
        return df, md
    return df


