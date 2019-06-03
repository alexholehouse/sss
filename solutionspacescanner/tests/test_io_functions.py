"""
Unit and regression test for io_functions.py
"""

# Import package, test suite, and other packages as needed
import solutionspacescanner
from solutionspacescanner import io_functions
import pytest
import sys

def test_solutionspacescanner_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "solutionspacescanner" in sys.modules


def test_readfile():
    """Sample test, will always pass so long as import statement worked"""

    fn = solutionspacescanner.get_data('look_and_say.dat')
    x = io_functions.readfile(fn)
    assert len(x) == 15
    
    
    

                          
