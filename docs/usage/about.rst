.. _About:

About
=========================================================


**SolutionSpaceScanner** is a software tool for generating customized ABSINTH parameter files which can be used to carry out solution space scanning using the CAMPARI simulation engine. Specifically, this is achieved using the command-line tool ``sss``, an approach that was recently applied to study the solution-dependence of intrinsically disordered regions by `Alex S. Holehouse <https://www.holehouselab.com/>`_ and `Shahar Sukenik <https://www.sukeniklab.com/>`_ [1].


.. _ABSINTH:

The ABSINTH implicit solvent model
***********************************

`CAMPARI <http://campari.sourceforge.net/>`_ is an all-atom simulation engine that facilitates both Monte Carlo and molecular dynamics simulations [2]. CAMPARI provides the reference implementation for the ABSINTH implicit solvent model, a model in which peptide-solvent interactions are described using a direct mean-field (DMFI) model [3]. The ABSINTH potential energy function (which we refer to as the Hamiltonian for convenience) describes the energy function that defines the instantaneous potential energy of a given configuration - put simply, one can take a three-dimensional configuration of a protein, run it "through" the ABSINTH Hamiltonian, and obtain a value in kcal/mol that reflects the potential energy of that particular conformation.

The ABSINTH Hamiltonian is defined as:

:math:`E_{total}=W_{solv}+U_{LJ}+W_{el}+U_{corr}`

Here the :math:`W_{solv}` term reflects the DMFI and provides the energetic contribution associated with the solvation of the set of chemical groups that exist in the polypeptide. The reference implementation of the ABSINTH model assumes the polypeptide is in an aqueous environment, and the associated energies associated with the DMFI are parameterized as such using experimentally measured values. However, these transfer free energies can be modified to reflect changes to the solution, allowing the 'chemical' nature of the underlying implicit solvent model to be changed in an arbitrary manner. 

The expression above  defines the general structure of the ABSINTH potential energy function. Each of these terms makes reference to a large number of parameters that specify the underlying molecular detail (ideal bond lengths, bond angles, partial charges, and group-wise reference free energies of solvation).  These parameters are provided by an ABSINTH parameter file. If that parameter file is modulated then the underlying interactions that occur in ABSINTH simulations can be altered. The SolutionSpaceScanner package was developed to automate the modulation of the solvation parameters within an ABSINTH parameter file.


Solution space scanning
*************************

The intrinsic tunability built into the ABSINTH model allows us to modulate the underlying free energies of solvation for different groups, allowing us to create chemically distinct solution conditions. This lets us ask how different types of solution environments may influence the conformational biases associated with a given polypeptide. For example, we can make a solution in which the backbone-solvent interactions become much less favorable than they would otherwise be (a more positive reference free energy of solvation).

We recently proposed that the systematic variation of the solution environment could be used as a way to fundamentally alter the conformational behavior of intrinsically disordered proteins. In much the same way that saturation mutagenesis allows for the role of sequence on conformational behaviour to be interrogated via a "sequence space scan", varying the solution environment allows for a "solution space scan". This approach was leveraged to examine how distinct solutions influence the behaviour of different intrinsically disordered proteins [3].


References
***********

[1] Holehouse, A.S., and Sukenik, S. (2019) Controlling Structural Bias in Intrinsically Disordered Proteins Using Solution Space Scanning (bioRxiv link: [to be posted])

[2] Vitalis, A., and Pappu, R.V. (2009). Methods for Monte Carlo Simulations of Biomacromolecules. In Annual Reports in Computational Chemistry, R.A. Wheeler, ed. (Elsevier), pp. 49–76.

[3] Vitalis, A., and Pappu, R.V. (2009). ABSINTH: A new continuum solvation model for simulations of polypeptides in aqueous solutions. J. Comput. Chem. 30, 673–699.


