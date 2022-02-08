"""
Utility functions shared by other modules
"""
import os
from pathlib import Path

import gisutils
import numpy as np
from scipy import interpolate
from shapely.geometry import MultiPolygon, Polygon
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
    print('grouping observations by geographic regions...')
    md = metadata.copy()
    if geographic_groups is not None:
        # get the group name from a dictionary
        geo_group_dict = {}
        if isinstance(geographic_groups, dict):
            for group_name, polygon in geographic_groups.items():
                # in the case of a shapefile, apply all shapes to that group name
                if isinstance(polygon, str) or isinstance(polygon, Path):
                    group_info = shp2df(str(polygon), dest_crs=metadata_crs)
                    # convert to polygons if needed
                    geoms = [Polygon(g) for g in group_info.geometry]
                    geo_group_dict[group_name] = MultiPolygon(geoms)
                else:
                    geo_group_dict[group_name] = polygon
        else:
            # option to supply shapefile with multiple groups/polygons
            if isinstance(geographic_groups, str) or isinstance(geographic_groups, Path):
                geographic_groups = [geographic_groups]
            for item in reversed(geographic_groups):
                group_info = shp2df(str(item), dest_crs=metadata_crs)
                groups = dict(zip(group_info[geographic_groups_col],
                                  group_info['geometry']))
                geo_group_dict.update(groups)
        for group_name, polygon in geo_group_dict.items():
            within = [g.intersects(polygon) for g in md.geometry]
            md.loc[within, 'geo_group'] = group_name
            print(f'{group_name}: {np.sum(within)} obs ({np.sum(within)/len(within):.1%})')
    return md


def cull_data_to_active_area(data, active_area, active_area_id_column=None,
                             active_area_feature_id=None,
                             data_crs=None, metadata=None):
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
        print('Culling {} sites outside of the model area defined by {}.'
              .format(np.sum(~within), active_area))
    df = df.loc[df_within]
    if metadata is not None:
        return df, md
    return df


