🔧 Installation
===============

🐍 PyPi
^^^^^^^

**EnzyMM** is hosted on GitHub, but the easiest way to install it is to download
the latest release from its `PyPi repository <https://pypi.python.org/pypi/pyjess>`_.
It will install all dependencies including the library of catalytic templates:

.. code:: console

    $ pip install --user enzymm

.. tip::

    To let **EnzyMM** read compressed files, write `.parquet`` tables
    and transformation matrices, install it including optional dependencies with:

    .. code::

        $ pip install --user enzymm[polars]


Conda
^^^^^

**EnzyMM** is also available as a `recipe <https://anaconda.org/bioconda/enzymm>`_
in the `bioconda <https://bioconda.github.io/>`_ channel. To install, simply
use the ``conda`` installer. Optionally, if you wish to search CIF structures,
install `gemmi <https://github.com/project-gemmi/gemmi/>`_ too:

.. code:: console

    $ conda install -c bioconda enzymm gemmi


.. Arch User Repository
.. ^^^^^^^^^^^^^^^^^^^^

.. A package recipe for Arch Linux can be found in the Arch User Repository
.. under the name `python-pyjess <https://aur.archlinux.org/packages/python-pyjess>`_.
.. It will always match the latest release from PyPI.

.. Steps to install on ArchLinux depend on your `AUR helper <https://wiki.archlinux.org/title/AUR_helpers>`_
.. (``yaourt``, ``aura``, ``yay``, etc.). For ``aura``, you'll need to run:

.. .. code:: console

..     $ aura -A python-pyjess


.. BioArchLinux
.. ^^^^^^^^^^^^

.. The `BioArchLinux <https://bioarchlinux.org>`_ project provides pre-compiled packages
.. based on the AUR recipe. Add the BioArchLinux package repository to ``/etc/pacman.conf``:

.. .. code:: ini

..     \[bioarchlinux\]
..     Server = https://repo.bioarchlinux.org/$arch

.. Then install the latest version of the package and its dependencies with ``pacman``:

.. .. code:: console

..     $ pacman -S python-pyjess


🐋 Docker container
^^^^^^^^^^^^^^^^^^^

You can also run **EnzyMM** via a `Docker <https://www.docker.com/>`_ container.
A container for every tagged release can be downloaded from the GitHub Container Repository.
Download the latest with:

.. code:: console

    $ docker pull ghcr.io/rayhackett/enzymm:latest


🖼️ Apptainer container
^^^^^^^^^^^^^^^^^^^^^^

You can also run **EnzyMM** via an `Apptainer <https://apptainer.org/>`_ container.
A container for every tagged release can be downloaded via ORAS from the GitHub Container Repository.
Download the latest with:

.. code:: console

    $ apptainer pull oras://ghcr.io/rayhackett/enzymm:latest

GitHub + ``pip``
^^^^^^^^^^^^^^^^

If, for any reason, you prefer to download the library from GitHub, you can clone
the repository and install the repository by running (with the admin rights):

.. code:: console

    $ git clone --recursive https://github.com/rayHackett/enzymm
    $ pip install --user ./enzymm

.. caution::

    Keep in mind this will install always try to install the latest commit,
    which may not even build, so consider using a versioned release instead.


GitHub + ``installer``
^^^^^^^^^^^^^^^^^^^^^^

If you do not want to use ``pip``, you can still clone the repository and
run ``build`` manually, although you will need to install the build 
dependencies.
dependencies are:
- rich
- PyJess (if you want to process cif files you will need PyJess[cif] or gemmi)
- readerwriterlock

.. code:: console

    $ git clone --recursive https://github.com/rayHackett/enzymm
    $ cd enzymm
    $ python -m build .

.. Danger::

    Installing packages without ``pip`` is strongly discouraged, as they can
    only be uninstalled manually, and may damage your system.