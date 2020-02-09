"""
Unit and regression test for io_functions.py
"""

# Import package, test suite, and other packages as needed
import solutionspacescanner
from solutionspacescanner import fos_calcs, io_functions
import pytest
import sys
import random



def test_get_explicit_phi():
    """
    Test we can read in a file - expect there to be 15 lines in this default file
    """

    #fos_calcs.get_explicit_phi(

    input_SGs = ['ALA','LYS']
    (SGs_to_modify, SG_count_dictionary) = io_functions.identify_used_residues('ALPLALPALPAPLPALPA', input_SGs)
    offset_dict = {'ALA':-0.2,
                   'LEU':0,
                   'PRO':0,
                   'PEP_PRO_BB':0,
                   'PEP_BB':0}

    print(SGs_to_modify)
    x = fos_calcs.get_explicit_phi(SG_count_dictionary, SGs_to_modify, offset_dict)
    
    assert x[0] == 1.9458810547025769
    assert x[1] == -61.668723126728466
    assert x[2] == -62.86872312672846
    assert len(x) == 3

        

    
