=====
Usage
=====

To use Double Decoupling Method (DDM) in a project

.. code-block:: console

    import ddm

In your terminal, use the command line

.. code-block:: console

    $ ddm --config your_config_file --complex pdb_complex


Example for a configuration file::

    [main]
    # Host code
    host=
    # Guest code
    guest=
    # Destination
    dest=

    [pick-reference]
    # Anchor points
    P1=
    P2=
    P3=
    L1=
    L2=
    L3=

This example contains the required parameters.

- Host code and Guest code are, respectively, the code of the protein and the code of the ligand in the pdb_complex file.
- Destination is the directory where outputs should be written.
- Anchor points are the atoms selected to monitors the movements of the ligand in respect of the receptor. See the section :ref:`pick_reference` for more information on how to pick them.

  If you don't specify those value, the program will stop when they are required. You can then look at the REFERENCE.pdb to pick them, add the values to your config file and relaunch the program using the same command.


More options can be added:

[main]:

- ff_parameters: Force field parameters for gromacs. Default is 'charmm36-jul2017.ff'.
- temperature: Set the temperature for the experiment. Default is '???'.
