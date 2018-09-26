=====
Usage
=====

To use Double Decoupling Method (DDM) in a project

.. code-block:: console

    from ddm.ddm import DDM
    ddm = DDM(config_file, pdb_complex)
    ddm.run()

If you want don't want to perform all the steps, you can use the individual modules in ddm.classes.
However, some steps required that others have already been performed (see :ref:`diagram`).

In your terminal, use the command line

.. code-block:: console

    $ ddm --config your_config_file --complex pdb_complex


Example for a configuration file::

    [main]
    # Host code
    host =
    # Guest code
    guest =
    # Destination
    dest =

    [pick-reference]
    # Anchor points
    P1 =
    P2 =
    P3 =
    L1 =
    L2 =
    L3 =

    [vba-unbound]
    # symmetry numbers to compute the symmetry correction
    symmetry_numbers = sigma_l, sigma_p, sigma_pl

This example contains the required parameters.

- Host code (**host**) and Guest code (**guest**) are, respectively, the code of the protein and the code of the ligand in the pdb_complex file.
- Destination (**dest**) is the directory where outputs should be written.
- Anchor points are the atoms selected to monitors the movements of the ligand in respect of the receptor. See the section :ref:`pick_reference` for more information on how to pick them.
    If you don't specify those value, the program will stop when they are required. You can then look at the REFERENCE.pdb to pick them, add the values to your config file and relaunch the program using the same command.
- **symmetry_numbers** for the ligand, the receptor and the complex to compute the correction due to the symmetry. If those numbers are not provided, the symmetry correction is not computed, then is equal to 0.


More options can be added:

[main]:

- **ff_parameters**: Force field parameters for gromacs. Default is charmm36-jul2017.ff.
- **temperature**: Set the temperature for the experiment. Default is ???.

[vba-bound]:

- **windows**: Factors to apply on constrain consecutively for dG vba-bound calculation. Default is 0.001, 0.01, 0.1, 0.2, 0.5, 1.0.
- **int_meth**: Integration method; WHAM or TI. Default is TI.
- **WHAM_PATH**: Required if WHAM has been chosen. The path for wham command. Default is ''.

