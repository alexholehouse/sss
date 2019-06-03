"""
Unit and regression test for io_functions.py
"""

# Import package, test suite, and other packages as needed
import solutionspacescanner
from solutionspacescanner import io_functions
import pytest
import sys
import random

def test_solutionspacescanner_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "solutionspacescanner" in sys.modules


def test_readfile():
    """
    Test we can read in a file - expect there to be 15 lines in this default file
    """

    fn = solutionspacescanner.get_data('look_and_say.dat')
    x = io_functions.readfile(fn)
    assert len(x) == 15


def test_error_exit():
    """
    Test error_exit exits and returns 1 as error code

    """
    with pytest.raises(SystemExit) as e:
        io_functions.error_exit('testing')
    assert e.type == SystemExit
    assert e.value.code == 1


def test_parse_residue_string():
    """
    Test residue_strinsg are parsed correctly

    """
    simple_AA = ['ALA','CYS','ASP','GLU','PHE','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']


    # first test empty string is appropiately dealt with (should error and exit)
    with pytest.raises(SystemExit) as e:
        io_functions.parse_residue_string('')
    assert e.type == SystemExit
    assert e.value.code == 1
    
    # next test invalid also trigger an exit with (1)
    with pytest.raises(SystemExit) as e:
        io_functions.parse_residue_string('ABA_DOG_CAT')
    assert e.type == SystemExit
    assert e.value.code == 1


    # next test invalid also trigger an exit with (1) because
    # we ignore glycine
    with pytest.raises(SystemExit) as e:
        io_functions.parse_residue_string('ABA_DOG_CAT_GLY')
    assert e.type == SystemExit
    assert e.value.code == 1

    # however should be able to extract ALA
    r = io_functions.parse_residue_string('ABA_DOG_CAT_ALA')
    assert len(r) == 1
    assert r[0] == 'ALA'

    # check all normal AAs
    for r in simple_AA:
        r_out = io_functions.parse_residue_string(r)
        assert r_out[0] == r

    # check peptide backbone
    r_out = io_functions.parse_residue_string('PEP-BB')
    assert r_out[0] == 'PEP_BB'

    # check HIS is correctly converted in into D/E/P versions
    r_out = io_functions.parse_residue_string('HIS')
    for i in ['HIE', 'HIP','HID']:        
        assert i in r_out


    for i in range(0,20):
        
        rset = random.choices(simple_AA, k=5)
        rset_string = "_".join(rset)
        r_out = io_functions.parse_residue_string(rset_string)
        for j in rset:
            assert j in r_out

    
    
        

        

                          
