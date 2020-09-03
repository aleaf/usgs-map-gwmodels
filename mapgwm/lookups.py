"""Module for importing lookup tables with information that covers the MERAS region
(stored in the lookups folder).

e.g. so that the contents of lookups/streamflow_site_comids.csv
can be import as a python dictionary with
>>> from mapgwm.lookups import streamflow_site_comids
"""
import os
import pandas as pd
import yaml


# get the path to this module
path, _ = os.path.split(__file__)
# location of lookup files relative to this module
inflows_lookup_file = os.path.join(path, 'lookups', 'streamflow_site_comids.csv')
nhdplus_vs_flowline_routing = os.path.join(path, 'lookups', 'nhdplus_v2_flowline_routing.csv')
aquifer_codes_data = os.path.join(path, 'lookups', 'aquifer_codes.yml')


def get_nhdplus_v2_flowline_routing():
    """ read in the NHDPLUS v2 flow routing table
        from a csv specified in the program
        /mapgwm/lookups/nhdplus_v2_flowline_routing.csv

        User adds::

            from mapgwm.lookups import flowline_routing

        to a script to retrieve flowline_routing as a dict
    """
    lookup = pd.read_csv(inflows_lookup_file)
    lookup = dict(zip(lookup.site_no, lookup.comid))
    return lookup


def get_streamflow_site_comids(group=None, groups=None):
    """ read in the site comid ID table 
        from a csv specified in the program
        /mapgwm/lookups/streamflow_site_comids.csv

        User adds::

            from mapgwm.lookups import streamflow_site_comids

        to a script to retrieve streamflow_site_comids as a dict
    """
    lookup = pd.read_csv(inflows_lookup_file)
    if groups is not None:
        lookup = lookup.loc[lookup.group.isin(groups)]
    if group is not None:
            lookup = lookup.loc[lookup.group == groups]
    lookup = dict(zip(lookup.site_no, lookup.comid))
    return lookup

def get_aq_cd_names():
    with open(aquifer_codes_data, 'r') as AQ:
        lookup = yaml.safe_load(AQ)
    return lookup


# execute the functions to read the data
# (each function only gets executed if its output variable is imported)
streamflow_site_comids = get_streamflow_site_comids()
flowline_routing = get_nhdplus_v2_flowline_routing()
aq_codes_dict = get_aq_cd_names()