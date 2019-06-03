#!/usr/bin/env python
"""
sss.py is the executable script that provides a command-line tool for performing 

"""

import os
import numpy as np
import argparse

# import local configuration files
from solutionspacescanner.configs import VALID, VALID_2_ONELETTER, SPACER_1
from solutionspacescanner import fos_calcs
from solutionspacescanner import io_functions
import solutionspacescanner


############################################################################################
if __name__=="__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-f", help="Input parameter file (requires an ABSINTH parameter file)") 
    parser.add_argument("-r", help="Residue string (format ALA_CYS_ASP_GLU etc.") 
    parser.add_argument("--fos_offset", help="Offset free energy of solvation being applied to each residue") 
    parser.add_argument("--fos_value",  help="Absolute value of the free energy of solvation being applied to each residue") 
    parser.add_argument("--fos_percentage",  help="Percentage delta WSM (between -100 and +100) value. If set MUST also provide a sequence")     
    parser.add_argument("--sequence",  help="Protein sequence used to calculate percentage delta WSM ")     
    parser.add_argument("-o", help="Output filename") 

    print("------------------------------------------------")
    print("  SolutionSpaceScanner (sss) version [%s] (%s)              " % (solutionspacescanner.__version__, solutionspacescanner.release_date))
    print("------------------------------------------------")

    # parse all passed arguments and set assumed default FOS_MODE
    args = parser.parse_args()
    FOS_MODE = 0

    # parse input file (ABSINTH parameter file)
    if (not args.f):
        print("Please pass an input ABSINTH parameter file ( -f <FILENAME> )")
        exit(1)
    else:
        # read the ABSINTH parameter file in - note no sanity checking is done when reading this file in... 
        # this is just to check we can read the file OK...
        content = io_functions.readfile(args.f)
        del content

    # parse residue string
    if not args.r:
        print("Please pass an input residue string (-r ALA_CYS_ASP_GLU etc.")
        print("Each residue with three letter code separated by an underscore")
        exit(1)
    else:
        
        # note - ALL sanity checking done inside parse_residue_string and will error
        # and exist automatically if things go awry
        residues = io_functions.parse_residue_string(args.r)
            

    ##
    ## rewire_fos lets you re-set the free-energy of solvation for one or more groups in three possibe ways.
    ## These are offset_mode, value_mode, and percentage_mode


    # --------------------------------------------------------------------------------------------------
    # If --fos_offset is used, we define a value that is applied as a fixed offset to each group defined in the residue
    # string
    if args.fos_offset:

        try:
            FOS_offset = float(args.fos_offset)
        except:
            io_functions.error_exit("Could not convert the --FOS_offset value to a number [%s] " % args.fos_offset)

        # finish up...
        FOS_value = 0 # set to default
        updated_residues = residues
        FOS_MODE = 1
        print("Using FOS offset mode. Each residue will have an offset of %5.5f applied" % (FOS_offset))


    # --------------------------------------------------------------------------------------------------
    # If --fos_values is used, we define an absolute value (in kcal/mol) that is set as the FOS value
    if args.fos_value:
        if FOS_MODE == 1:
            io_functions.error_exit("Passed both a FOS value and a FOS offset - please use only one")

        # finish up...
        FOS_value = float(args.fos_value)
        FOS_offset = 0 # set to default
        FOS_MODE = 2
        updated_residues = residues

    # --------------------------------------------------------------------------------------------------
    # If --fos_percentage is used we calculate an fixed offset value that is determined such that the new
    # MTFE is changed from the 'aqeous' MTFE by this percentage. To do this we need to know what the protein
    # sequences is, as the MTFE is calculated based on the amino acids present
    if args.fos_percentage:
        if FOS_MODE == 1 or FOS_MODE == 2:
            io_functions.error_exit("Passed both a FOS percentage and a FOS offset or value - please use only one")

        # finish up with FOS_MODE=3
        if not args.sequence:
            io_functions.error_exit("Must provide a valid amino acid sequence when using --fos_percentage")
            
        FOS_MODE = 3

        updated_residues = fos_calcs.identify_used_residues(args.sequence, residues)
        FOS_offset = fos_calcs.compute_value_from_percentage(args.sequence, float(args.fos_percentage), updated_residues)        
        FOS_value = 0 # default



    # if FOS_MODE was never set (i.e. no --fos_* was used)
    if FOS_MODE == 0:
        io_functions.error_exit("Please pass either a FOS (free energy of solvation) offset (--fos_offset) or a FOS value (--fos_value)")

    # if we passed an outpfile set this, else use a default value for the output 'rewired' parameter file
    if args.o:
        outname = args.o
    else:
        outname = 'new_params.prm'

    # write out a new ABSINTH parameter file
    io_functions.update_ABSINTH_parameter_file(args.f, outname, updated_residues, FOS_MODE, FOS_value, FOS_offset)


    
    
        
