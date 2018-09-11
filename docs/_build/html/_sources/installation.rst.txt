.. highlight:: shell

============
Installation
============

From sources
------------

The sources for Double Decoupling Method (DDM) can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/donadef/ddm


Once you have a copy of the source (it is recommanded to create a virtual environment, see _`Handle virtual environmement` for more infos), you can install it with:

.. code-block:: console

    $ pip isntall -r requirements.txt
    $ python setup.py install

If you want to be able to launch calculation using only the 'ddm' command, use:

.. code-block:: console

    $ pip install --editable .


.. _Github repo: https://github.com/donadef/ddm


Handle virtual environment
--------------------------

To create a virtual environment for your project, use `pyenv`_ (python version management) and `virtualenv`_ (manage virtual environments).

.. _pyenv: https://github.com/pyenv/pyenv
.. _virtualenv: https://github.com/pyenv/pyenv-virtualenv


You can follow the installation instruction from the GitHub repositories. Here’s a resume if you have an Ubuntu distribution:

pyenv:
++++++

+ Clone the repository in a directory ~/.pyenv (as recommanded):

.. code-block:: console

    $ git clone https://github.com/pyenv/pyenv.git ~/.pyenv


+ Define environment variables:

.. code-block:: console

    $ echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
    $ echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc

+ Add pyenv init to your shell

.. code-block:: console

    $ echo -e 'if command -v pyenv 1>/dev/null 2>&1; then\n  eval "$(pyenv init -)"\nfi' >> ~/.bashrc

+ Install python versions (ddm is written in python 3)

.. code-block:: console

    $ pyenv install 3.6.5

If you encounter dependency problems, check pyenv's `wiki`_.

.. _wiki: https://github.com/pyenv/pyenv/wiki/common-build-problems

virtualenv:
+++++++++++

+ Clone the repository into pyenv's plugin directory:

.. code-block:: console

    $ git clone https://github.com/pyenv/pyenv-virtualenv.git $(pyenv root)/plugins/pyenv-virtualenv

+ Create a virtual environment for ddm:

.. code-block:: console

    $ pyenv virtualenv 3.6.5 venv_ddm

+ Activate the virtual environment:

.. code-block:: console

    $ pyenv activate venv_ddm

The prompt will indicate that the environment is active.

+ Deactivate the environment:

.. code-block:: console

    $ pyenv deactivate

+ Delete the virtual environment:

.. code-block:: console

    $ pyenv uninstall venv_ddm


