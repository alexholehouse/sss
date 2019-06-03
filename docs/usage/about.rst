About
=========================================================
**SolutionSpaceScanner** is a software tool for generating customized ABSINTH parameter files which can be used to carry out solution space scanning (sss) using the CAMPARI simulation engine. 

SolutionSpaceScanner implements a method used in the paper "*Controlling Structural Bias in Intrinsically Disordered Proteins Using Solution Space Scanning*" by `Alex S. Holehouse <https://www.holehouselab.com/>`_ and `Shahar Sukenik <https://www.sukeniklab.com/>`_.

CAMPARI is an all-atom simulation engine that facilitates both Monte Carlo and molecular dynamics simulations [1]. CAMPARI provides the reference implementation for the ABSINTH implicit solvent paradigm, a model in which peptide-solvent interactions are described using a direct mean-field model (DMFI). The ABSINTH Hamiltonian describes the energy function that defines the instantaneous potential energy of a given configuration. 

The ABSINTH Hamiltonian is defined as:

:math:`E_{total}=W_{solv}+U_{LJ}+W_{el}+U_{corr}`

Here 

To do a developmental install, type

``pip install -e .``

Dependencies
************************

You need to install ``numpy``


Subheader
--------------------------
Sample subheader



References
--------------------------
[1] Vitalis, A., and Pappu, R.V. (2009). Methods for Monte Carlo Simulations of Biomacromolecules. In Annual Reports in Computational Chemistry, R.A. Wheeler, ed. (Elsevier), pp. 49–76.

[2] Vitalis, A., and Pappu, R.V. (2009). ABSINTH: A new continuum solvation model for simulations of polypeptides in aqueous solutions. J. Comput. Chem. 30, 673–699.
