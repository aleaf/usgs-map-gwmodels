
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

import warnings
from packaging import version
import gisutils

if version.parse(gisutils.__version__) < version.parse('0.2.4'):
    warnings.warn('USGS-MAP-gwmodels requires gis-utils >= 0.2.4'
                  '\nPlease pip install --upgrade gis-utils')