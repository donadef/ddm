.. highlight:: shell

============
Installation
============


Stable release
--------------

To install Double Decoupling Method (DDM), run this command in your terminal:

.. code-block:: console

    $ pip install ddm

This is the preferred method to install Double Decoupling Method (DDM), as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for Double Decoupling Method (DDM) can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/donadef/ddm

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/donadef/ddm/tarball/master

Once you have a copy of the source (it is recommanded to create a virtual environment), you can install it with:

.. code-block:: console

    $ pip isntall -r requirements.txt
    $ python setup.py install

If you want to be able to launch calculation using only the 'ddm' command, use:

.. code-block:: console

    $ pip install --editable .


.. _Github repo: https://github.com/donadef/ddm
.. _tarball: https://github.com/donadef/ddm/tarball/master
