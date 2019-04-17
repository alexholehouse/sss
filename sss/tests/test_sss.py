"""
Unit and regression test for the sss package.
"""

# Import package, test suite, and other packages as needed
import sss
import pytest
import sys

def test_sss_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "sss" in sys.modules
