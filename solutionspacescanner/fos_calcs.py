from . import configs 
from . import offset_calculator
from . import io_functions


# ******************************************************************
def compute_value_from_percentage(target_percentage, SG_count_dictionary, updated_SGs_to_modify):
    """
    Function that, considering the full amino acid sequence ($sequence), the desired change
    in MTFE percentage ($percentage) and the set of residues that can be changed ($residues)
    determines the offset that should be evenly applied to each residue type defined in $residues
    to achieve this percentage change. This is achived by calling a function in offset_calculator
    package    

    """

    # call the 
    offset = offset_calculator.get_delta_percentage_W_solvmax(target_percentage, SG_count_dictionary, updated_SGs_to_modify)

    return offset


def get_explicit_phi(SG_count_dictionary, updated_SGs_to_modify, offset_dict):
    """

    """

    wt_wsolvmax      =  offset_calculator.get_offset_W_solv_max(SG_count_dictionary, updated_SGs_to_modify, 0)
    altered_wsolvmax =  offset_calculator.get_offset_W_solv_max(SG_count_dictionary, updated_SGs_to_modify, offset_dict)    
    phi = 100*( (altered_wsolvmax - wt_wsolvmax)/ wt_wsolvmax)
    
    return [phi, wt_wsolvmax, altered_wsolvmax]
    
        


def bonus_phi(SGs_to_modify, seq, FOS_offset):
    """
    SGs_to_modify is a list of solvation groups
    
    seq is an amino acid sequence
    
    FOS_offset is a dictionary of SGs to FOS offset values, where each offset will be applied 

    """
    
    (updated_SGs_to_modify, SG_count_dictionary) = io_functions.identify_used_residues(seq, SGs_to_modify)

    # compute the delta % W_{solv}^{max} based on
    output = get_explicit_phi(SG_count_dictionary, updated_SGs_to_modify, FOS_offset)

    print('+-------------------------------------------+')
    print('|     Summary of MTFE on W_solv^max         |')
    print('+-------------------------------------------+')
    print('')
    print("sss is changing the W_solv_max from %3.3f to %3.3f" % (output[1], output[2]))
    print("The associated percentage change in phi is:")
    print("\nDelta Phi: %1.3f %%\n" %(output[0]))
    print("This change is ONLY true given the the passed sequence:\n%s"%(seq))
    print('-------------------------------------------')
    print('')
    

