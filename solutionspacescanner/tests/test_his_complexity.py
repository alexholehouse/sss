import solutionspacescanner
from solutionspacescanner import io_functions, configs
import pytest
import sys
import random


def test_parse_histidine_correctly():
    """
    Test residue_strings are parsed correctly

    """

    # make sure HIE/HID/HIP in MTFE file raise exceptions
    for i in ['HIE','HID','HIP']:
        fn = solutionspacescanner.get_data('test_data/1MUREA_%s.mtfe'%(i))
        with pytest.raises(io_functions.SSSException):
            FOS_offset = io_functions.parse_mtfe_file(fn)
