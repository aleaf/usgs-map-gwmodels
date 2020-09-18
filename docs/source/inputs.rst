========================================
Additional description of input data
========================================

Head Observation Input
-----------------------
Description of input fields for groundwater level times series:

    .. literalinclude:: ../../mapgwm/tests/data/headobs/GW_monthly_stats_test.txt
        :start-after: # -----
        :end-before: # -----

Description of input fields for groundwater level metadata:

    .. literalinclude:: ../../mapgwm/tests/data/headobs/GW_monthly_meta_test.txt
        :start-after: # -----
        :end-before: # -----

Aquifer Code Names
-----------------------
Local aquifer codes are mapped to the following names

    .. literalinclude:: ../../mapgwm/lookups/aquifer_codes.yml
        :language: yaml
        :start-after: aquifer_code_names:
        :end-before: # mapping
        :dedent: 2

Regional Aquifer Code Names
----------------------------------
For the purposes of the MAP project, local aquifer codes
are generalized to the following names that represent the equivalent unit(s)
across the Mississippi Embayment:

    .. literalinclude:: ../../mapgwm/lookups/aquifer_codes.yml
        :language: yaml
        :start-after: regional_aquifer:
        :dedent: 2
