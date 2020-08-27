"""Module for importing lookup tables with information that covers the MERAS region
(stored in the lookups folder).

e.g. so that the contents of lookups/streamflow_site_comids.csv
can be import as a python dictionary with
>>> from mapgwm.lookups import streamflow_site_comids
"""
import os
import pandas as pd


# get the path to this module
path, _ = os.path.split(__file__)
# location of lookup files relative to this module
inflows_lookup_file = os.path.join(path, 'lookups', 'streamflow_site_comids.csv')
nhdplus_vs_flowline_routing = os.path.join(path, 'lookups', 'nhdplus_v2_flowline_routing.csv')


def get_nhdplus_v2_flowline_routing():
    lookup = pd.read_csv(inflows_lookup_file)
    lookup = dict(zip(lookup.site_no, lookup.comid))
    return lookup


def get_streamflow_site_comids(group=None, groups=None):
    lookup = pd.read_csv(inflows_lookup_file)
    if groups is not None:
        lookup = lookup.loc[lookup.group.isin(groups)]
    if group is not None:
            lookup = lookup.loc[lookup.group == groups]
    lookup = dict(zip(lookup.site_no, lookup.comid))
    return lookup


# execute the functions to read the data
# (each function only gets executed if its output variable is imported)
streamflow_site_comids = get_streamflow_site_comids()
flowline_routing = get_flowline_routing()