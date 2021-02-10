## sss - a package for performing computational solution space scanning
##
## Written by Alex Holehouse 
## Developed by Alex Holehouse and Shahar Sukenik
## See LICSENCE for copyright information
##
## offset_calculator.py
## Contains functionality 
##
##

import numpy as np
import os
import errno
import argparse
from . import configs


################################################################################################
####                                                                                        ####
####                      Global parameters                                                 ####
####                                                                                        ####
################################################################################################


## ...........................................................................
##
def get_offset_W_solv_max(SG_count_dictionary, updated_SGs_to_modify, offset):
    """


    Function that given an dictionary defining the number of each type of SG in the 
    system  (${SG_count_dictionar}), the set of SGs to be updated (updated_SGs_to_modify)
    and an offset that can be either a dictionary of values or a single fixed value, this
    function calculates the W_solv_max for the protein associated with the SG_count_dictionary.

    Note there is some logic involved here, which includes

    (1) Correcting backbone SASA based on the presence of a sidechain (note we don't right now
        correct the sidechain SASA for the presence of a backbone, but assume is not a major contributor
        - this may be something to address later

    (2) the peptide backbone is dealt with explicitly by counting the number of proline and non-proline
        residues and then calculating the backbone contribution accordingly. Because glycine is not explicitly
        included we also calculate the number of glycines based on TOTAL number of backbone moeities minus
        number of sidechain-continaing SGs.

    Parameters
    ...........
    SG_count_dictionary : dict
        Dictionary where key-value pairs are solvation group (SG) and count.

    updated_SGs_to_modify : list
        Non-redundant list of SGs that are a subset of the SGs in SG_count_dictionary.
    
    offset : float OR dict
        Either a fixed value applied to every SG in in updated_SGs_to_modify list OR a dictionary
        where EVERY SG in SG_count_dictionary MUST be able to look up what's going on!
    

    """
    max_fos = 0.0
    non_gly = 0


    # for each solvation group in the protein
    for SG in SG_count_dictionary:

        # allows this function to take offset as EITHER
        # a fixed value or a fixed value
        if isinstance(offset, dict):
            local_offset = offset[SG]
        else:
            local_offset = offset

        # .......................................................................................
        # FIRST compute sidechain contribution
            
        # we deal with non-proline backbones explicitly below, while proline backbone is dealt with 
        if SG == 'PEP_BB' or SG == 'PEP_PRO_BB':
            continue

        # if that group is in the set to modify (note that this includes proline sidechain here)
        if SG in updated_SGs_to_modify:            
            max_fos = max_fos + SG_count_dictionary[SG]*(configs.FOS_baseline[SG] + local_offset);
        else:
            max_fos = max_fos + SG_count_dictionary[SG]*(configs.FOS_baseline[SG])

        # .......................................................................................
        # backbone for non-glycines
            
        # compute a correction factor - this is the SASA of the backbone of the SG of interest
        # divided by the SASA of a GLY backbone residue
        BB_correction_factor = (configs.SASA_MAX[SG][1]/configs.SASA_MAX['GLY'][1]) # BB correction for presence of sidechain
        
        # PROLINE
        # if we found a proline... rest of this is about proline backbone stuff
        if SG == 'PRO':

            # update offset IF a variable offset was provided (if not stick with fixed value)
            if isinstance(offset, dict):
                local_offset = offset['PEP_PRO_BB']

            # if we are meant to updated PEP_PRO_BB...
            if 'PEP_PRO_BB' in updated_SGs_to_modify:
                max_fos = max_fos + SG_count_dictionary['PEP_PRO_BB']*(configs.FOS_baseline['PEP_PRO_BB'] + local_offset)*BB_correction_factor
            else:
                max_fos = max_fos + SG_count_dictionary['PEP_PRO_BB']*(configs.FOS_baseline['PEP_PRO_BB'])*BB_correction_factor

        # NON PROLINE
        # else we found a non-proline residue which must have a backbone residue
        else:

            non_gly = non_gly + SG_count_dictionary[SG]

            # update offset IF a variable offset was provided
            if isinstance(offset, dict):
                local_offset = offset['PEP_BB']

            # if we are meant to updated PEP_BB - note, both of these compuet the addition to max_fos based on
            # (1) number of residues we're on (SG), multiplied by the FOS_baseline value with/without offset
            # multiplied by the BB_correction factor
            if 'PEP_BB' in updated_SGs_to_modify:
                max_fos = max_fos + SG_count_dictionary[SG]*(configs.FOS_baseline['PEP_BB'] + local_offset)*BB_correction_factor
            else:
                max_fos = max_fos + SG_count_dictionary[SG]*(configs.FOS_baseline['PEP_BB'])*BB_correction_factor

    # .......................................................................................
    # glycine            
    # finally we deal with glycine - note we have the 'IF' statement here BASICALLY as the edge
    if 'PEP_BB' in SG_count_dictionary: 
        ngly = SG_count_dictionary['PEP_BB'] - non_gly

        if isinstance(offset, dict):
            local_offset = offset['PEP_BB']
        else:
            local_offset = offset

        # note no BB_correction_factor because glycine BB are fully exposed by definition
        if 'PEP_BB' in updated_SGs_to_modify:
            max_fos = max_fos + ngly*(configs.FOS_baseline['PEP_BB']+local_offset)
        else:
            max_fos = max_fos + ngly*(configs.FOS_baseline['PEP_BB'])

    return max_fos
        




## ...........................................................................
##
def get_delta_percentage_W_solvmax(target_percentage, SG_count_dictionary, updated_SGs_to_modify):
    """

    Parameters
    ----------------

    target_percentage : float
        Target percentage we wish the percentage Delta W_{solv}^{max} to shift by

    SG_count_dictionary: dictionary
        Already parsed dictionary that contains key:value pairs of solvation group and count for 
        the sequence

    updated_SGs_to_modify : list
        List that is a subset of the solvation groups in the SG_count_dictionary



    """

    if target_percentage == 0:
        print("No offset required (target percentage change = 0)")
        return 0

    print("Input percentage: %3.3f " % target_percentage)


    # fist compute the W_{solv}^{max} for the sequence under aqeous 
    # conditions
    wt_W_solv_max = get_offset_W_solv_max(SG_count_dictionary, updated_SGs_to_modify, 0)

    W_solv_max_resolution = 0.005

    trajectory = []
    offset = -W_solv_max_resolution


    if target_percentage < 0:
        PC = -1
    else:
        PC = 1
    target_percentage=abs(target_percentage)

    delta_percentage_W_solv_max = 100*( (get_offset_W_solv_max(SG_count_dictionary, updated_SGs_to_modify, offset) - wt_W_solv_max)/wt_W_solv_max)

    while delta_percentage_W_solv_max < target_percentage:
        offset = offset - W_solv_max_resolution
        delta_percentage_W_solv_max = 100*( (get_offset_W_solv_max(SG_count_dictionary, updated_SGs_to_modify, offset) - wt_W_solv_max)/wt_W_solv_max)
        trajectory.append(delta_percentage_W_solv_max)

    # check which of the two values that straddle the boundary is closed and choose the best of the two
    
    if abs(trajectory[-2]-target_percentage) < abs(delta_percentage_W_solv_max - target_percentage):

        # the + here means IF we find that the penultimate delta_percentage_W_solv_max is better we 
        # step back and make the offset 
        offset = offset + W_solv_max_resolution
        delta_percentage_W_solv_max = trajectory[-2]


    # if we have a negative we are reducing the RFOS 
    if PC < 0:
        offset = -offset
        delta_percentage_MTFE = -delta_percentage_W_solv_max

    print("Offset of %4.3f to groups %s gives a percent of %1.3f (original = %5.1f kcal/mol, rewired = %5.1f kcal/mol) " %(offset, updated_SGs_to_modify, delta_percentage_W_solv_max, wt_W_solv_max, get_offset_W_solv_max(SG_count_dictionary, updated_SGs_to_modify, offset)))


    return offset






## ...........................................................................
##
def get_mtfe_derived_W_solv_max(FOS_offset, original_fos, solvation_group_count, PRO_PEP):
    """
    FOSS_offset is the off

    """
    original_fos = 0.0

    solvation_groups = list(original_fos.keys())

    BB_count = 0
    BB_PRO_count = 0

"""
    # for each sidechain
    for sg in solvation_groups:
        
        #if sg == ''

        c = solvation_group_count[sg]

        original_fos = original_fos + c*original_fos[sg]

        if sg == 'PRO' and c >0:
            BB_PRO_count = c
        else:
            BB_count = BB_count + c


    original_fos = original_fos
            

original_[AA_IDX]*(FOS_baseline[INT2AA[AA_IDX]]+offset);
        
        AA_IDX=AA2INT[AA]
        
        # if sidechain group is being changed
        if AA in group:
            pass

            # sidechain max possible FOS with offset set
            
        else:
            pass

            # sidechain max possible FOS (no offset)
            max_fos = max_fos + resvector[AA_IDX]*FOS_baseline[INT2AA[AA_IDX]];

        # if backbone not going to be changed then deal with this here
        if 'B' not in group:

            if PRO_PEP and AA == 'P':
                BB_correction_factor = 1 # BB correction for presence of proline - assume that no correction needed
                max_fos = max_fos + resvector[AA_IDX]*pro_backbone_baseline*BB_correction_factor
                
            else:                
                # and backbone max possible FOS (no offset) with fractional correction cos sidechains take up space!
                BB_correction_factor = (SASA_MAX[AA_IDX][1]/SASA_MAX[5][1]) # BB correction for presence of sidechain
                max_fos = max_fos + resvector[AA_IDX]*backbone_baseline*BB_correction_factor

    # if backbone is being corrected then we DID NOT include it in the previous loop and we deal woth it here..
    if 'B' in group:

        for AA in AAs:        
            AA_IDX=AA2INT[AA]
        
            # sidechains have already been dealt with so if 'B' was in group we are now JUST adjusting the backbone

            # and backbone max possible FOS (no offset)
            if PRO_PEP and AA == 'P':
                BB_correction_factor = 1 # BB correction for presence of proline - assume that no correction needed
                max_fos = max_fos + resvector[AA_IDX]*(pro_backbone_baseline+offset)*BB_correction_factor
            else:
                BB_correction_factor = (SASA_MAX[AA_IDX][1]/SASA_MAX[5][1]) # BB correction for presence of sidechain
                max_fos = max_fos + resvector[AA_IDX]*(backbone_baseline+offset)*BB_correction_factor
    
    return max_fos
"""
