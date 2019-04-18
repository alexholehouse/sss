from . configs import VALID, VALID_2_ONELETTER, SPACER_1
from sss import offset_calculator


# ******************************************************************
def compute_value_from_percentage(sequence, percentage, residues):
    """
    Function that, considering the full amino acid sequence ($sequence), the desired change
    in MTFE percentage ($percentage) and the set of residues that can be changed ($residues)
    determines the offset that should be evenly applied to each residue type defined in $residues
    to achieve this percentage change. This is achived by calling a function in offset_calculator
    package    

    """

    # cyle through each valid residue, convering into the relevant one-letter
    # code to create a groupstring, which is fed into the offset_calculator
    # function
    groupstring=""
    for r in residues:
        groupstring = groupstring+VALID_2_ONELETTER[r]
        
        
    # print some status information...
    print("")
    print("Full sequence         : %s" % sequence)
    print("Residues to be changed: %s" % groupstring)


    # check to ensure at least one of the residues we want to change is found
    # in the sequence
    no_residues=True

    # edgecase for backbone
    if 'B' in groupstring:
        no_residues=False
        
    # extract backbone from checks as iff B is in the groupstring
    # it is obvously present in the sequence!
    groupstring_minus_bb=groupstring.replace('B','')

    notfound=''
    for i in groupstring_minus_bb:
        if i in sequence:
            no_residues=False
        else:
            notfound=notfound+i

    if no_residues:
        error_exit("Using FOS percentage mode but none of the passed residues [%s] are found in the sequence [ %s ] "% (", ".join(list(groupstring)), sequence))
    else:
        if len(notfound)>0:
            print("[WARNING]: Using FOS percentage mode to change several residue groups, but the following are not found in this sequence: %s" % (",".join(list(notfound))))
            print('           ONLY the residues found in the sequence will be changed')
            print("")
            

    # call the get_deleta_percentage_MTFE function
    offset = offset_calculator.get_delta_percentage_MTFEs(sequence, percentage, groupstring)

    return offset

# ******************************************************************
def identify_used_residues(sequence, residues):
    """
    Function that, considering the full amino acid sequence ($sequence), the desired change
    in MTFE percentage ($percentage) and the set of residues that can be changed ($residues)
    determines the offset that should be evenly applied to each residue type defined in $residues
    to achieve this percentage change. This is achived by calling a function in offset_calculator
    package    

    """

    # cyle through each valid residue, convering into the relevant one-letter
    # code to create a groupstring, which is fed into the offset_calculator
    # function

    # for each residue passed
    updated_residues=[]
    firstime=True
    for r in residues:        
        single_letter = VALID_2_ONELETTER[r]
        if single_letter != 'B':
            if single_letter in sequence.upper():
                updated_residues.append(r)
            else:
                if firstime:
                    print('\n')
                    firstime=False
                print('***** [WARNING]: Using percentage MTFE but residue %s not found in sequence - removing from altered residue set...' % (r))
        else:
            updated_residues.append(r)
            
    # if firsttime is false means one or more values were written, so add a trailing newline
    if firstime is False:
        print('\n')
            
    return updated_residues

        
