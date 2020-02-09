## sss - a package for performing computational solution space scanning
##
## Written by Alex Holehouse 
## Developed by Alex Holehouse and Shahar Sukenik
## See LICSENCE for copyright information
##
## configs.py
## configs contains information used by function within the ss package. 
## This includes validation info, standard ABSINTH free energy of solvation
##

AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

# spacer for writing ABSINTH parameter files
SPACER_1 = '              '

THREELETTER_TO_SOLAVTION_GROUP = {'ALA':'ALA',
                                  'CYS':'CYS',
                                  'ASP':'ASP',
                                  'GLU':'GLU',
                                  'PHE':'PHE', 
                                  'HIS':'HIE',
                                  'HIE':'HIE',
                                  'HID':'HID',
                                  'HIP':'HIP',
                                  'ILE':'ILE',
                                  'LYS':'LYS',
                                  'LEU':'LEU',
                                  'MET':'MET',
                                  'ASN':'ASN',
                                  'PRO':'PRO',
                                  'GLN':'GLN',
                                  'ARG':'ARG',
                                  'SER':'SER',
                                  'THR':'THR',
                                  'VAL':'VAL',
                                  'TRP':'TRP',
                                  'TYR':'TYR',
                                  'PEP-BB':'PEP_BB',
                                  'PEP-PRO-BB':'PEP_PRO_BB'}

ONELETTER_TO_SOLVATION_GROUP  = {'A':'ALA',
                                'C':'CYS',
                                'D':'ASP',
                                'E':'GLU',
                                'F':'PHE', 
                                'H':'HIE',
                                'I':'ILE',
                                'K':'LYS',
                                'L':'LEU',
                                'M':'MET',
                                'N':'ASN',
                                'P':'PRO',
                                'Q':'GLN',
                                'R':'ARG',
                                'S':'SER',
                                'T':'THR',
                                'V':'VAL',
                                'W':'TRP',
                                'Y':'TYR'}

# baseline sidechain FOS values as taken from the ABINSTH implicit solvent model. Note
# the FOS of G is set to 0 here (no sidechain). Values here in kcal/mol.
FOS_baseline = {'ALA': 1.9,
                'CYS': -1.2,
                'ASP': -107.3,
                'GLU': -107.3,
                'PHE': -0.8,
                'HIP': -95.0,
                'HIE': -10.3,
                'HID': -10.3,
                'ILE': 2.2,
                'LYS': -100.9,
                'LEU': 2.3,
                'MET': -1.4,
                'ASN': -9.7,
                'PRO': 2.0,
                'GLN': -9.4,
                'ARG': -100.9,
                'SER': -5.1,
                'THR': -5.0,
                'VAL': 2.0,
                'TRP': -5.9,
                'TYR': -6.1,
                'PEP_BB': -10.1,
                'PEP_PRO_BB': -8.5}                  

# MAX SASA for sidechain and backbone in each residue, as measured from 
# ACE-XX-NME dipeptide (note HIE/HID/HIP all were HIE)
SASA_MAX = {}
SASA_MAX['ALA'] = [7.581871795654296875e+01, 7.607605743408203125e+01]   # A ALA
SASA_MAX['CYS'] = [1.154064483642578125e+02, 6.787722015380859375e+01]   # C CYS
SASA_MAX['ASP'] = [1.302558288574218750e+02, 7.182710266113281250e+01]   # D ASP
SASA_MAX['GLU'] = [1.617985687255859375e+02, 6.805746459960937500e+01]   # E GLU
SASA_MAX['PHE'] = [2.093871002197265625e+02, 6.598278808593750000e+01]   # F PHE
SASA_MAX['GLY'] = [0.000000000000000000e+00, 1.149752731323242188e+02]   # G GLY
SASA_MAX['HIE'] = [1.808149414062500000e+02, 6.750666809082031250e+01]   # H HIS
SASA_MAX['HID'] = [1.808149414062500000e+02, 6.750666809082031250e+01]   # H HIS
SASA_MAX['HIP'] = [1.808149414062500000e+02, 6.750666809082031250e+01]   # H HIS
SASA_MAX['ILE'] = [1.727196502685546875e+02, 6.034464645385742188e+01]   # I ILE
SASA_MAX['LYS'] = [2.058575897216796875e+02, 6.871156311035156250e+01]   # K LYS
SASA_MAX['LEU'] = [1.720360412597656250e+02, 6.451246643066406250e+01]   # L LEU
SASA_MAX['MET'] = [1.847660064697265625e+02, 6.778076934814453125e+01]   # M MET
SASA_MAX['ASN'] = [1.427441253662109375e+02, 6.680493164062500000e+01]   # N ASN
SASA_MAX['PRO'] = [1.342914733886718750e+02, 5.583909606933593750e+01]   # P PRO
SASA_MAX['GLN'] = [1.733262939453125000e+02, 6.660184478759765625e+01]   # Q GLN
SASA_MAX['ARG'] = [2.364875640869140625e+02, 6.673487854003906250e+01]   # R ARG
SASA_MAX['SER'] = [9.587133026123046875e+01, 7.287202453613281250e+01]   # S SER
SASA_MAX['THR'] = [1.309214324951171875e+02, 6.421310424804687500e+01]   # T THR
SASA_MAX['VAL'] = [1.431178131103515625e+02, 6.172962188720703125e+01]   # V VAL
SASA_MAX['TRP'] = [2.545694122314453125e+02, 6.430991363525390625e+01]   # W TRP
SASA_MAX['TYR'] = [2.225183105468750000e+02, 7.186695098876953125e+01]   # Y TYR

def listsetsort(x):
    x = list(x)
    x = set(x)
    x = list(x)
    x.sort()
    return x

def get_valid_solgroups():
    return listsetsort(THREELETTER_TO_SOLAVTION_GROUP.values())

def get_valid_three_letter_residues():
    return listsetsort(THREELETTER_TO_SOLAVTION_GROUP.keys())

def get_valid_one_letter_residues():
    return listsetsort(ONELETTER_TO_SOLAVTION_GROUP.keys())
