Quickstart
=========================================================


Installation
*************

Installation of the **solutionspacescanner** scanner package provides the commandlinetool ``sss``, a general tool for rewiring ABINSTH parameter files.

**solutionspacescanner** is available from PyPI. To install **solutionspacescanner** simply use `pip`:

	$ pip install solutionspacescanner

To ensure installation was succesfull you can run the following command from your terminal::

	$ sss --version

Usage
*************

Running the command ::

	$ sss --help

Should print the following information to the screen ::

	------------------------------------------------
	  SolutionSpaceScanner (sss) version [0.0.0+24.g7edb330.dirty] (November 2019)
	------------------------------------------------
	usage: sss [-h] [-f F] [-r R] [--fos_offset FOS_OFFSET]
    	       [--fos_value FOS_VALUE] [--fos_percentage FOS_PERCENTAGE]
        	   [--sequence SEQUENCE] [-o O] [--proline-pep]
	           [--mtfe_file MTFE_FILE] [--version]
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -f F                  Input parameter file (requires an ABSINTH parameter
	                        file)
	  -r R                  Residue string (format ALA_CYS_ASP_GLU etc.
	  --fos_offset FOS_OFFSET
	                        Offset free energy of solvation being applied to each
	                        residue
	  --fos_value FOS_VALUE
	                        Absolute value of the free energy of solvation being
	                        applied to each residue
	  --fos_percentage FOS_PERCENTAGE
	                        A percentage as defined as phi (between -100 and +100)
	                        value. If set MUST also provide a sequence
	  --sequence SEQUENCE   Protein sequence used to calculate phi
	                        
	  -o O                  Output filename
	  --proline-pep         If included, proline backbone is altered like the
	                        peptide backbone (default off)
	  --mtfe_file MTFE_FILE
	                        MTFE file that contains a mapping of residue to MTFE
	                        and a scaling factor
      --version             Print the current version string	  


Once installed ``sss`` allows for the re-wiring of an existing ABSINTH parameter files in a controlled and customizable way. Specifically, ``sss`` allows you to alter the reference free-energy of solvation in four distinct ways:

1. *Absolute mode* An absolute value for the reference free energy of solvation value associated with one or more amino acid groups is defined. This uses the ``--fos_value`` flag. FOS values are in kcal/mol. Note that every residue defined in the residue string will be set to the same value.

2. *Offset mode* An absolute offset will be applied to the reference free energy of solvation. Again this value is in kcal/mol.

3. *Percentage mode* The reference free energy of solvation for each solvation group is modulated by an absolutely equivalent amount such that the new OVERALL contribution from the maximally solvated protein (the maximum transfer free energy (defined as :math:`W^{max}_{solv}` in [1]) has changed by a fixed percentage value relative to the original value under aqeous conditions. Practically, this is done by first computing :math:`W^{max}_{solv}`  for the molecule under aqeous conditions, then calculating by how much :math:`W^{max}_{solv}` would need to change (in kcal/mol) to achieve the designated percentage change, and then distributing that difference out evenly across the designated solvation groups. This means two different solutions can be - from the perspective of the whole molecule - equivalent, but _how_ changes to the solvation groups are distributed will be very different. This is the mode used for all changes in [1].

4. Using an MTFE file. An MTFE file defines a fixed offset to be applied to every solvation group. For more information on the MTFE file format see below.


For a much more in-depth description of these three modes and how the **percentage mode** changes are calculated please see the supplementary information provided in [1].

The various command-line options are described in detail below.

``-h / --help`` : **help**, Prints the help information (as shown above)


``-f``  : **file**, defines the input ABSINTH parameter file that will be used as the base for the changes. Note that ``sss`` does not ensure a passed ABSINTH parameter file is a valid parameter file, although will ensure that all residues to be changed can be changed. Parameter files can be obtained with the full CAMPARI release (examples include ``abs3.2_opls.prm``).

``-R``  : **residues** defines the residue string for the set of residues that are to be altered. Specifically, a residue string is of the format ALA_CYS_ASP *etc.*, -i.e. three letter amino acid symbols separated by an underscore. The order these are provided is irrelevant, but they must be valid three-letter codes for one of the standard twenty amino acids. **Note** reference free energy of solvation for glycine cannot be changed without changing ALL of the backbone values. To avoid an inadvertent change to all backbone groups, if the user wishes to change the backbone the non-standard residue name PEP-BB must be provided (e.g. ALA_PEP-BB_CYS).



``--fos_value`` : **fos absolute value** - defines the absolute value to be used as the reference free energy of solvation for each of the residues defined in the residue string (by ``-r``)

``--fos_offset`` : **fos offset value** - absolute offset (can be positive or negative) to be applied to the reference free energies of solvation associated with the residues defined by ``-r``.

``--fos_percentage`` : **fos percentage value** - defines the percentage that the overall maximum transfer free energy of the entire molecule (:math:`W^{max}_{solv}`) will change by. 


Examples
*************

Make the hydrophobic residues more hydrophobic ::

	sss -f abs3.2_opls.prm --fos_offset 1.0 -r ALA_LEU_VAL_ILE_MET


Note here the file ``abs3.2_opls.prm`` is taken from the ABSINTH parameter files and should reflect an actual file that is present in the current directory. Assuming this works correctly this should print the following information to the screen::


	------------------------------------------------
	  SolutionSpaceScanner (sss) version [0.0.0+24.g7edb330.dirty] (November 2019)
	------------------------------------------------
	Using FOS offset mode. Each residue will have an offset of 1.00000 applied
	Updating residue ALA from 1.9 to  2.90
	Updating residue VAL from 2.0 to  3.00
	Updating residue LEU from 2.3 to  3.30
	Updating residue ILE from 2.2 to  3.20
	Updating residue MET from -1.4 to -0.40

	File [new_params.prm] written succesfully

As you can see, all the aliphatic hydrophobes defined in the input residue string have been shifted by 1 kcal/mol up in their free energy of solvation.

MTFE File format
****************

An MTFE file is a tab-separated input file for which each every solvation group should have a key-value pair. Comments can be included and should be defined using the '#' symbol. An aexample of a valid MTFE file is included below::

	#
	# Units in in cal/mol/res
	# Note scalar here is 5
	# These are comments that are ignored when the file is processed
	
	SCALAR  5
	ALA  0
	CYS  0
	ASP  0
	GLU  0
	PHE  -200   # you can also write comments inline
	HIS  -100
	ILE  0
	LYS  0
	LEU  0
	MET  0
	ASN  0
	PRO  0
	GLN  0
	ARG  0
	SER  0
	THR  0
	VAL  0
	TRP  -200
	TYR  -200
	PEP_BB  -50
	PEP_PRO_BB  -30

Note that EVERY single one of these solvation groups MUST be defined. Finally, the SCALAR keyword defines a fixed multiplier that is used when each of these values. The scalar is useful when calibrating experimental data with ABSINTH.

References
***********

[1] Holehouse, A.S., and Sukenik, S. (2020) Controlling Structural Bias in Intrinsically Disordered Proteins Using Solution Space Scanning, JCTC - online access (`Access article <https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.9b00604>`_)


.. [*] Note that ``solutionspacescanner`` will be added to PyPI in the near future.


