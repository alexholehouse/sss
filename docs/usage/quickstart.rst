Quickstart
=========================================================


Installation
*************

Installation of the **SolutionStateScanner** scanner package provides the commandlinetool ``sss``, a general tool for rewiring ABINSTH parameter files.

To install **SolutionStateScanner** download the early release candidate `zip file from here <https://www.holehouselab.com/s/solutionspacescanner.zip>`_ [*]_. Once downloaded, this file can be used to install **SolutionSpaceScanner** installed using ``pip`` by running ::

	$ pip install solutionspacescanner.zip

To ensure installation was succesfull you can run the following command from your terminal::

	$ sss --help 
	
This should print the following information to the screen ::

	------------------------------------------------
	  SolutionSpaceScanner (sss) version [<SOME VERSION>] (June 2019)
	------------------------------------------------
	usage: sss [-h] [-f F] [-r R] [--fos_offset FOS_OFFSET]
	           [--fos_value FOS_VALUE] [--fos_percentage FOS_PERCENTAGE]
	           [--sequence SEQUENCE] [-o O]
	
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
	                        Percentage delta WSM (between -100 and +100) value. If
	                        set MUST also provide a sequence
	  --sequence SEQUENCE   Protein sequence used to calculate percentage delta
	                        WSM
	  -o O                  Output filename
	  

Usage
************************

Once installed ``sss`` allows for the re-wiring of an existing ABSINTH parameter files in a controlled and customizable way. Specifically, ``sss`` allows you to alter the reference free-energy of solvation in three distinct ways:

1. *Absolute mode* An absolute value for the reference free energy of solvation value associated with one or more amino acid groups is defined. This uses the ``--fos_value`` flag. FOS values are in kcal/mol. Note that every residue defined in the residue string will be set to the same value.

2. *Offset mode* An absolute offset will be applied to the reference free energy of solvation. Again this value is in kcal/mol.

3. *Percentage mode* The reference free energy of solvation for each solvation group is modulated by an absolutely equivalent amount such that the new OVERALL contribution from the maximally solvated protein (the maximum transfer free energy (defined as :math:`W^{max}_{solv}` in [1]) has changed by a fixed percentage value relative to the original value under aqeous conditions. Practically, this is done by first computing :math:`W^{max}_{solv}`  for the molecule under aqeous conditions, then calculating by how much :math:`W^{max}_{solv}` would need to change (in kcal/mol) to achieve the designated percentage change, and then distributing that difference out evenly across the designated solvation groups. This means two different solutions can be - from the perspective of the whole molecule - equivalent, but _how_ changes to the solvation groups are distributed will be very different. This is the mode used for all changes in [1].

For a much more in-depth description of these three modes and how the **percentage mode** changes are calculated please see the supplementary information provided in [1].

The various command-line options are described in detail below.

``-h / --help`` : **help**, Prints the help information (as shown above)


``-f``  : **file**, defines the input ABSINTH parameter file that will be used as the base for the changes. Note that ``sss`` does not ensure a passed ABSINTH parameter file is a valid parameter file, although will ensure that all residues to be changed can be changed. Parameter files can be obtained with the full CAMPARI release (examples include ``abs3.2_opls.prm``).

``-R``  : **residues** defines the residue string for the set of residues that are to be altered. Specifically, a residue string is of the format ALA_CYS_ASP *etc.*, -i.e. three letter amino acid symbols separated by an underscore. The order these are provided is irrelevant, but they must be valid three-letter codes for one of the standard twenty amino acids. **Note** reference free energy of solvation for glycine cannot be changed without changing ALL of the backbone values. To avoid an inadvertent change to all backbone groups, if the user wishes to change the backbone the non-standard residue name PEP-BB must be provided (e.g. ALA_PEP-BB_CYS).



``--fos_value`` : **fos absolute value** - defines the absolute value to be used as the reference free energy of solvation for each of the residues defined in the residue string (by ``-r``)

``--fos_offset`` : **fos offset value** - absolute offset (can be positive or negative) to be applied to the reference free energies of solvation associated with the residues defined by ``-r``.

``--fos_percentage`` : **fos percentage value** - defines the percentage that the overall maximum transfer free energy of the entire molecule (:math:`W^{max}_{solv}`) will change by. 


References
***********

[1] Holehouse, A.S., and Sukenik, S. (2019) Controlling Structural Bias in Intrinsically Disordered Proteins Using Solution Space Scanning (bioRxiv link: [to be posted])


.. [*] Note that ``solutionspacescanner`` will be added to PyPI in the near future.


