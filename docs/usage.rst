=====
Usage
=====

To use Double Decoupling Method (DDM) in a project::

    import ddm

In your terminal, use the command line

.. code-block:: console

    $ ddm_cli --config your_config_file --complex pdb_complex


Example for a configuration file::

    [main]
    # Host code
    host=CB7
    # Guest code
    guest=AD0
    # Where to write the outputs ?
    dest=/storage/donadef/work/CB7-AD0

    [pick-reference]
    # Anchor points
    P1=115
    P2=120
    P3=124
    L1=136
    L2=128
    L3=129

This example contains the required options. More options can be added:

TODO
