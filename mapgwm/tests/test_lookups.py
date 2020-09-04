import pytest
from mapgwm.lookups import streamflow_site_comids
from mapgwm.lookups import aq_codes_dict
from mapgwm.lookups import flowline_routing

def test_lookups():
    assert(len(aq_codes_dict['aquifer_code_names']) > 0)
    assert(len(aq_codes_dict['regional_aquifer']) > 0)
    assert(isinstance(streamflow_site_comids, dict))
    assert(isinstance(flowline_routing, dict))


