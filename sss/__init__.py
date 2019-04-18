"""
SSS
Python package to perform solution state scanning and generate ABSINTH parameter files
"""

# Add imports here
from .sss import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


