
import solutionspacescanner
import pytest
import sys

test_data_dir = solutionspacescanner.get_data('test_data')

@pytest.fixture(scope='session', autouse=True)
def input_abs3p2(request):    
    """
    Function that returns a string that 

    """

    return '%s/abs3.2_opls.prm'%(test_data_dir)


