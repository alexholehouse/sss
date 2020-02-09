"""
Unit and regression test for io_functions.py
"""

# Import package, test suite, and other packages as needed
import solutionspacescanner
from solutionspacescanner import io_functions, configs
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


def test_readfile_abs(input_abs3p2):
    """
    Test we can read in a file - expect there to be 15 lines in this default file
    """

    fn = solutionspacescanner.get_data(input_abs3p2)
    x = io_functions.readfile(fn)

    # 3402 is how many lines in the abs3.2 file!
    assert len(x) == 3402


def test_error_exit():
    """
    Test error_exit exits and returns 1 as error code

    """
    with pytest.raises(SystemExit) as e:
        io_functions.error_exit('testing')
    assert e.type == SystemExit
    assert e.value.code == 1

def test_sanitize_sequence():
    """
    """

    # check bad ones fail
    for d in ['APLSBLPLAPLAL', 'LKLPLPLALALAP.', 'ASDA-ASDASD']:
        with pytest.raises(SystemExit) as e:
            io_functions.sanitize_sequence(d)
        assert e.type == SystemExit
        assert e.value.code == 1

    # check good ones work
    for d in ['asLPASLPAlal', 'AcDeFGhIKlMnPqRsTvWy']:
        out = io_functions.sanitize_sequence(d)
        assert out == d.upper()

def test_parse_mtfe_file():

    test_data_dir = solutionspacescanner.get_data('test_data')

    
    filenames = ['bad_mtfe.mtfe',
                 'bad_mtfe1.mtfe',
                 'bad_mtfe2.mtfe',
                 'bad_mtfe3.mtfe',
                 'bad_mtfe4.mtfe']

    # check that a bunch of examples with different errors fail
    for fn in filenames:
        full = '%s/%s'%(test_data_dir, fn)

        with pytest.raises(SystemExit) as e:
            io_functions.parse_mtfe_file(full)
        assert e.type == SystemExit
        assert e.value.code == 1
        
    # check a good one works too!
    full = '%s/%s'%(test_data_dir, 'good_mtfe.mtfe')
    x = io_functions.parse_mtfe_file(full)

    # expected length 
    assert len(x) == 23


def test_identify_used_residues():
    

    input_SGs = ['ALA','LYS']
    x = io_functions.identify_used_residues('ALPLALPALPAPLPALPA', input_SGs)
    assert len(x[0]) == 1
    assert x[0][0] == 'ALA'
    print(x[1])
    assert x[1]['ALA'] == 6
    assert x[1]['LEU'] == 6
    assert x[1]['PRO'] == 6
    assert x[1]['PEP_BB'] == 12
    assert x[1]['PEP_PRO_BB'] == 6


    input_SGs = ['PEP_BB','LYS']
    x = io_functions.identify_used_residues('GGGGGGGGGGGG', input_SGs)
    assert len(x[0]) == 1
    assert x[0][0] == 'PEP_BB'
    assert x[1]['PEP_BB'] == 12


    # can't distingsuih different types of HIS in sequence
    with pytest.raises(SystemExit) as e:
        input_SGs= ['HIS','HIE','HIP']
        x = io_functions.identify_used_residues('ALPLALPAHLPHAPLPALPA', input_SGs)
    assert e.type == SystemExit
    assert e.value.code == 1

    with pytest.raises(SystemExit) as e:
        input_SGs = ['VAL']
        x = io_functions.identify_used_residues('ALPLALPALPAPLPALPA', input_SGs)
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
    for i in ['HIE']:
        assert i in r_out


    for i in range(0,20):
        
        rset = random.choices(simple_AA, k=5)
        rset_string = "_".join(rset)
        r_out = io_functions.parse_residue_string(rset_string)
        for j in rset:
            assert j in r_out

    # check that the 'all' keyword works
    r_out = io_functions.parse_residue_string('all')
    for i in simple_AA:
        assert i in r_out
        
    assert 'PEP_BB' in r_out
    assert 'PEP_PRO_BB' in r_out
        
    
    
        

        

                          
