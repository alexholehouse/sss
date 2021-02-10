"""
io_functions

Code that deals with parsing/reading/reporting input and output


"""
from time import gmtime, strftime

from . import configs
from .sssexception import SSSException
from . import solutionspacescanner
import os


## ******************************************************************
##
##
def readfile(filename):
    """
    Function that takes a filestring as input and converts it into a list of lines,
    returning the list. Nothing fancy.

    """
    with open(filename) as fh:
        content = fh.readlines()
    return content



## ******************************************************************
##
##
def error_exit(msg):
    """

    Function that takes an error message, prints it in a nicely formatted way
    and then triggers an exit with a custom exception.


    """
    print("")
    print("###############################################################")
    print("\nERROR:  %s"  % str((msg)))
    print("")
    print("###############################################################")
    raise(SSSException('\n\nExiting gracefully. See error message above...'))



## ******************************************************************
##
##
def sanitize_sequence(s):
    """
    Function that parses protein string (one-letter amino acid sequence) and verifies each residue 
    can be correctly dealt with

    """

    seq = s.upper()
    for a in seq:
        if a not in configs.AAs:
            error_exit('---- Found invalid amino acid (%s)'%(a))
    return seq



## ******************************************************************
##
##
def parse_residue_string(rstring, PRO_PEP=False):

    """
    Function that takes an input residue string (ALA_CYS_ASP_...) and converts it into a list of valid
    solvation groups. 

    Notable functionality:

        * GLY is ignored (for the peptide backbone to be changed the input value PEP-BB must be provided)
        * All possible histidine codes (HIS/HIE/HID/HIP) are treated in the same way and universally changed 
              HIE/HID/HIP (the three types of HIS in ABSINTH
        * Only valid amino acid codes are processed - non-valid ones are ignored with a warning
        * Valid amino acid codes are defined in configs.py file within the solutionspacescanner package


    Parameters
    ..........

    rstring : string
        The rstring is a string that seperates three-letter coded amino acids by an underscore, and is the 
        raw input string provided to the sss tool (-r flag).

    PRO_PEP : bool
        Flag that if true means we treat PEP_PRO_BB like PEP_BB.


    Returns
    .......
    
    list
        Returns a list of valid solvation that we are going to try and change. The groups that are 
        actually changed will be a subset of these and is limited to the intersection of those SGs defined
        in this function and those present in the 
    
    
    """    

    # firstly we extract the valid set of solvation groups and the valid set of three-letter
    # residues from .configs. Note these are both non-redundant ordered lists 
    valid_SGs              = configs.get_valid_solgroups()
    valid_3letter_residues = configs.get_valid_three_letter_residues()

    # if 'all' passed we're using ALL solvation groups. Note the remove lines
    # here mean that the valid_SGs contain only HIE as the one 
    if rstring == 'all':
        print("Using residue grouping: 'all'\n")
        valid_SGs.remove('HIP')
        valid_SGs.remove('HID')
        return valid_SGs
        
    # split the residue string up using '_'
    split_string = rstring.strip().upper().split("_")

    # define the list of SGs we will return 
    return_SGs = []

    # for each residue name as extracted from the XXX_XXX_XXX string...
    for r in split_string:

        # if glycine is passed warn we can't use glycine and continue
        if r == 'GLY':
            print("[WARNING]: GLYCINE cannot be changed without changing ALL the backbone. To avoid unintented consequences, please use PEP-BB instead of GLY. This residue is being ignored")
            continue

        # if residue is not in the valid set of three-letter residue identifiers
        elif r not in valid_3letter_residues:
            print("[WARNING]: Residue [%s] is not recognized as a valid residue" %(r))
            continue

        # if residue is any of the 'his' type residues add all three HIS-related solvation
        # groups (note HIS is not a solvation group!)
        elif r == 'HIS':
            print("[INFO]:    Defaulting HIS->HIE")
            # note campari doesn't use HIS, so neutral histadine is generally represented as HIE 
            return_SGs.append('HIE')
            continue

        # else lookup in the conversion table
        else:
            return_SGs.append(configs.THREELETTER_TO_SOLAVTION_GROUP[r])

    # check we found at least one residue...
    if len(return_SGs) == 0:
        error_exit("No viable residues were provided to change [%s] ..."%(rstring))

    # if this falg was true, and IF we had added the peptide backbone then ALSO add the
    # proline backbone
    if PRO_PEP:
        if 'PEP_BB' in return_SGs:            
            return_SGs.append('PEP_PRO_BB')

    # finally remove any duplicates and sort for consistent behaviour
    return_SGs = set(return_SGs)
    return_SGs = list(return_SGs)
    return_SGs.sort()
    
    return return_SGs



## ******************************************************************
##
##
def identify_used_residues(sequence, SGs_to_modify):
    """
    Takes in an amino acid sequence and a set of groups to modify and 

    (1) updates the set of groups to ensure ONLY groups that are present in
        the sequence are considered and

    (2) calculates a dictionary mapping the number of each SG to its official
        SG name


    For some further explanation:

    This code is ONLY run when computing a percentage Phi value. The SGs_to_modify
    is a list of parsed solvation groups as built from the residue string 
    (i.e. -r ALA_CYS_HIS...). This string is actually parsed by the function
    parse_residue_string().

    The seqeuence variable is the one-letter amino acid string for the amino acid
    sequence.

    This function filters out the SGs_to_modify based on the 

    

    Parameters
    ----------

    sequence : string
        Amino acid string.

    SGs_to_modify : list
        This is a list of official solvation groups (SGs) that sss has been
        passed via the -r flag to modify (parsed from 3-letter codes into
        bona fide solvation groups already)

    Returns
    -------

    tuple : (updated_SGs_to_modify, SG_count_dictionary)

        updated_SGs_to_modify - this is an updated list that's a union or subset
        of the input list SGs_to_modify, but has removed ANY SGs that are not present
        in the residue

        SG_count_dictionary - this is a dictionary of ALL solvation groups in the protein,
        essentially containing all the sequence information decomposed into a more convenient
        format for further work.   

    """

    # cyle through each valid residue, convering into the relevant one-letter
    # code to create a groupstring, which is fed into the offset_calculator
    # function

    # local function
    # ..............................
    def addtodict(d,k):
        if k not in d:
            d[k] = 0
        d[k] = d[k]+1
        return d
    # ..............................

    ##
    ## STAGE 1: Count ALL solvation groups within the polypeptide
    ##
    SG_count_dictionary={}
    for position in sequence:
        
        # if gly incrememnt backbone and continue
        if position == 'G':
            SG_count_dictionary = addtodict(SG_count_dictionary, 'PEP_BB')
            continue

        # if proline incrememnt proline backbone and proline sidechain - note this means we 
        # do NOT increment the PEP_BB when we get a proline
        elif position == 'P':
            SG_count_dictionary = addtodict(SG_count_dictionary, 'PRO')
            SG_count_dictionary = addtodict(SG_count_dictionary, 'PEP_PRO_BB')            
            continue

        elif position == 'H':

            # if multiple types of histidine were passed to modify then this should fail.
            # This bit of code counts if any of the histidine three-letter codes were passed (note that HIS->HIE in the parse_residue_string() function
            # so the ONLY possible histidine codes here will be HIE, HID, or HIP
            local = 0
            for d in ['HIE','HID','HIP']:
                if d in SGs_to_modify:
                    local = local + 1
            if local > 1:
                error_exit('Currently there is no way to distinguish between different histidine residues, please only try and modify one of HIS/HIE/HID/HIP')

            # if exactly one type of histidine was passed we're going to assume ALL histdines in the sequence are of this type and modify all
            elif len(set(SGs_to_modify).intersection(['HIE','HID','HIP']))  == 1:

                # get the SG
                sg = list(set(SGs_to_modify).intersection(['HIE','HID','HIP']))[0]
                print('INFO: Assuming all HIS residues in the sequence are of the solvation group: %s' % (sg))
                SG_count_dictionary = addtodict(SG_count_dictionary, sg)
                SG_count_dictionary = addtodict(SG_count_dictionary, 'PEP_BB')        
                
            # else a HIS is in the sequence but we're not modifying histidine, so let's just use the default based on the
            # OENLETTER_TO_SOLVATION_GROUP mapping (i.e. sg will be HIE)
            else:
                pass
                    
        # convert amino acid residue into a pair of solvation groups 
        # (BB and sidechian) and add these
        sg = configs.ONELETTER_TO_SOLVATION_GROUP[position]    
        SG_count_dictionary = addtodict(SG_count_dictionary, sg)
        SG_count_dictionary = addtodict(SG_count_dictionary, 'PEP_BB')        

    ##
    ## STAGE 2: Cross-reference 
    ##
    #print(SGs_to_modify)
    updated_SGs_to_modify = []
    firstime = True # this is actually for nice warning formatting...
    for sg in SGs_to_modify:
        if sg in SG_count_dictionary:
            updated_SGs_to_modify.append(sg)
        else:
            if firstime:
                print('\n')
                firstime=False
            print('***** [WARNING]: Solvation group %s was not found in sequence - removing this group from the set of residues to alter for the purposes of calculating Phi'%(sg))
                                                    
    # if firsttime is false means one or more values were written, so add a trailing newline
    if firstime is False:
        print('\n')

    # IF none of the modified residues are left, throw and error!
    if len(updated_SGs_to_modify) == 0:
        error_exit("Using FOS percentage mode but none of the passed residues are found in the sequence")
                    
    return (updated_SGs_to_modify, SG_count_dictionary)



## ******************************************************************
##
##
def _discard_comment(line):
    """ 
    Function that discards anything after the first comment symbol, which
    for ABSINTH parameter files is a '#'.

    For example, line 

    "fos ALA -10 # this is a comment"

    would return
    
    "fos ALA -10"
    

    Parameters
    ..........

    line : string
        Input line from an ABSINTH parameter file

    Returns
    .......

    Same string with comments removed. 

    """
    return line.strip().split('#')[0]



## ******************************************************************
##
##    
def get_original_fos_values(filename):
    # for now this is hardcoded!!
    return configs.FOS_baseline



## ******************************************************************
##
##    
def parse_mtfe_file(mtfe_file):
    """
    Function that takes in an MTFE file and returns a dictionary that maps three-letter codes back to MTFE modification and a scalar that acts as a uniform multiplier for 
    each value when apploed to the system.

    Some implementation details:
        Duplicate values in the input file trigger an exception and exit

        If name - value paris are found where the name is not in the list of allowed values the line is ignored

        If value cannot be cast to a float this triggers an exception and exit

        If a name is missing the function will fail 
        

    Parameters
    -------------
    mtfe_file : string
        String that defines the location of an MTFE file, a text file with the following format

        # comments start with a '#' symbol
         
        <name>  <value>

        Where name-value paris reflect amino acid name and a signed value applied to those residues.
    
        Nineteen residues excluding GLY must be defined, as must be PEP_BB, PEP_PRO_BB and SCALAR as three additional names

    Returns
    -------------
    dictionary 
    
         Returns a dictionary of names to values that is guarenteed to be 23 elements in length

    

    """

    # EVERY one of these must be explicitly defined
    required_names = ['SCALAR',  'ALA',  'CYS',  'ASP',  'GLU',  'PHE',  'HIS',   'ILE',  'LYS',  'LEU',   'MET',   'ASN',   'PRO',   'GLN',   'ARG',   'SER',   'THR',   'VAL',   'TRP',   'TYR',   'PEP_BB',  'PEP_PRO_BB']


    MTFE_RETURN={}
    for i in required_names:
        MTFE_RETURN[i] = False
    
    content = readfile(mtfe_file)
    

    print('-------------------------------------------------------')
    print('Parsing MTFE input file [%s]'%(mtfe_file))
    for line in content:

        # if a fully comment line, skip
        if len(line.strip()) < 1:
            continue

        if line.strip()[0] == '#':
            continue
        
        c = _discard_comment(line)
        
        c_split = c.split()
        if len(c_split) == 2:
            
            if c_split[0] in required_names:
                if MTFE_RETURN[c_split[0]] is not False:
                    error_exit('Found duplicate of %s in MTFE input file. Exiting to be safe...'%(c_split[0]))

                try:
                    MTFE_RETURN[c_split[0]] = float(c_split[1])
                except Exception:
                    error_exit('Unable to convert value to float on line %s'%(line))
            else:
                print('%s not in expected inputs... skipping'%(c_split[0]))

    for i in MTFE_RETURN:
        if MTFE_RETURN[i] is False:
            error_exit('Input MTFE file did not define %s'%(i))

    S = MTFE_RETURN['SCALAR']
    final_mtfe_values={}

    for i in MTFE_RETURN:
        if i == 'SCALAR':
            print('SCALAR factor set to %3.3f' % (MTFE_RETURN[i]))
            continue
        else:
            # convert cal/mol/res to kcal/mol/res
            final_mtfe_values[i] = (MTFE_RETURN[i]*S)/1000.0
            print ('Offset %s %3.3f' %(i, final_mtfe_values[i]))

    # FINALLY we add in additional offsets for derived residues
    final_mtfe_values['HIE'] = final_mtfe_values['HIS'] 
    final_mtfe_values['HID'] = final_mtfe_values['HIS'] 
    final_mtfe_values['HIP'] = final_mtfe_values['HIS'] 
    
    # remove HIS (not used in parameter file)
    final_mtfe_values.pop('HIS')

    return final_mtfe_values



## ******************************************************************
##
##
def update_ABSINTH_parameter_file(inname, outname, residues, FOS_MODE, FOS_value, FOS_offset):
    """
    Function that will update an ABSINTH parameter file with the new FOS values as defined in FOS_value or
    FOS_offset depending on FOS_MODE.

    Note this function does some sanity checking to make sure everything looks good!

    """
    
    # define these for sanity checking once we're finished with the file..
    sanity_check = {}
    for i in residues:
        sanity_check[i] = 0
        

    # read in parameter file
    content = readfile(inname)

    # open a new file and then cycle through each line of the previously read
    with open(outname, 'w') as fh:
        fh.write("### This file was generated by sss version %s  on  %s \n \n" % (solutionspacescanner.__version__, strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        if FOS_MODE == 4:
            fh.write("### Specifically, using input MTFE file\n")
            

        # for each line in the original ABSINTH parameter file...
        for line in content:

            # remove the comment and split the line by whitespace
            l = _discard_comment(line)
            sl = l.split()

            # if this line contains parameter information...
            if len(sl) == 3:

                # if the file is a fos line...
                if sl[0].lower() == 'fos':

                    groupname = sl[1].upper().strip()

                    # if the residue in question is one of the ones we want to change...                                        
                    if groupname in residues:
                        
                        if sanity_check[groupname] > 0:
                            error_exit('\n\nError: Found two examples of residues [%s] in the parameter file [%s]' % (sl[1].upper(), inname))
                            
                        # update the sanity-check counter. We use this to ensure we don't find multiple copies of the same protein
                        sanity_check[groupname] = sanity_check[groupname] + 1
                            
                        # change the value using the approriate approach depending on mode, and the
                        # print that this has changed
                        if FOS_MODE == 2:
                            new_val = FOS_value
                        elif FOS_MODE == 1 or FOS_MODE == 3:
                            new_val = float(sl[2].upper()) + FOS_offset
                        elif FOS_MODE == 4:
                            new_val = float(sl[2].upper()) + FOS_offset[groupname]
                            
                        print("Updating residue %s from %s to %5.2f" % (sl[1].upper(), sl[2], new_val ))
                        fh.write("fos%s%s %5.2f # UPDATED by sss \n" %(configs.SPACER_1, sl[1].upper(), new_val))

                    else:
                        fh.write(line)
                else:
                    fh.write(line)
            else:
                fh.write(line)
    
    # .. end of with clause ...
    
    for i in sanity_check:
        if sanity_check[i] == 0:

            # delete the previously generated file to avoid any possible hint that it worked,
            # print out the error and exit
            os.remove(outname)
            error_exit('\n\nDid not find residue [%s] in the parameter file, exiting...' %(i))

    print("")        
    print("File [%s] written succesfully" % outname)
    print("")        
    print("")        

