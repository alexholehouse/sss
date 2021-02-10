"""
solutionspacescanner.py
Python package to perform solution state scanning and generate ABSINTH parameter files

Handles the primary functions
"""

release_date = 'February 2021'
from . _version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    pass
