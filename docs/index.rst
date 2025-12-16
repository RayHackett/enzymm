.. matcher documentation master file, created by
   sphinx-quickstart on Thu Aug 22 13:20:10 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



.. .. |AUR| image:: https://img.shields.io/aur/version/python-pyjess?logo=archlinux&style=flat&maxAge=3600
..    :target: https://aur.archlinux.org/packages/python-pyjess
..    :class: dark-light

EnzyMM - The Enzyme Motif Miner |Stars|
=======================================

.. |Stars| image:: https://img.shields.io/github/stars/rayhackett/enzymm.svg?style=social&label=Star&maxAge=3600
   :target: https://github.com/rayhackett/enzymm/stargazers
   :class: dark-light

*Identify catalytic sites in protein structures with the Enzyme Motif Miner!*

|Actions| |Coverage| |PyPI| |Wheel| |Docker| |Apptainer| |Versions| |versions| |License| |Source| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/RayHackett/enzymm/test.yml?branch=main&style=flat&maxAge=300
   :target: https://github.com/RayHackett/Enzymm/actions/workflows/test.yml
   :class: dark-light

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/rayhackett/enzymm?logo=codecov&style=flat&maxAge=3600
   :target: https://codecov.io/gh/rayhackett/enzymm/
   :class: dark-light

.. |PyPI| image:: https://img.shields.io/pypi/v/enzymm.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/enzymm
   :class: dark-light

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/enzymm?style=flat&maxAge=3600
   :target: https://anaconda.org/bioconda/enzymm
   :class: dark-light

.. |Wheel| image:: https://img.shields.io/pypi/wheel/enzymm?style=flat&maxAge=3600
   :target: https://pypi.org/project/enzymm/#files
   :class: dark-light

.. |Docker| image:: https://img.shields.io/badge/Docker-GHCR-blue?logo=docker
   :target: https://github.com/users/rayhackett/packages/container/package/enzymm
   :class: dark-light

.. |Apptainer| image:: https://img.shields.io/badge/Apptainer-SIF-blue?logo=apptainer&style=flat
   :target: https://github.com/rayhackett/enzymm/releases/latest
   :class: dark-light

.. |Versions| image:: https://img.shields.io/pypi/pyversions/enzymm.svg?style=flat&maxAge=3600
   :target: https://pypi.org/project/enzymm/#files
   :class: dark-light

.. |versions| image:: https://img.shields.io/github/v/tag/rayhackett/enzymm?label=version&sort=semver
   :target: https://github.com/rayhackett/enzymm/tags
   :class: dark-light

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/
   :class: dark-light

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=3600&style=flat
   :target: https://github.com/rayhackett/enzymm/
   :class: dark-light

.. |Issues| image:: https://img.shields.io/github/issues/rayhackett/enzymm.svg?style=flat&maxAge=600
   :target: https://github.com/rayhackett/enzymm/issues
   :class: dark-light

.. |Docs| image:: https://img.shields.io/readthedocs/enzymm?style=flat&maxAge=3600
   :target: http://enzymm.readthedocs.io/en/stable/?badge=stable
   :class: dark-light

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-green.svg?maxAge=3600&style=flat
   :target: https://github.com/rayhackett/enzymm/blob/main/CHANGELOG.md
   :class: dark-light

.. |Downloads| image:: https://static.pepy.tech/personalized-badge/enzymm?period=total&units=INTERNATIONAL_SYSTEM&left_color=grey&right_color=GREEN&left_text=downloads
   :target: https://pepy.tech/projects/enzymm
   :class: dark-light

Overview
--------

**EnzyMM** is a Python tool to identify known patterns of catalytic residues in protein structures.
.. We also provide a webserver where you can run **EnzyMM** on single queries
.. `here <https://www.ebi.ac.uk/thornton-srv/m-csa/enzymm>`_

Catalytic sites, which define the function of proteins, are among the most conserved
parts of enzyme structures. It makes sense therefore, that their similarities can be
identified even if the rest of the protein has diverged throughout the course of
evolutionary history. Similarly doing so immediately gives biological insight into
possible functions.
Proteins carry out a myriad of different functions. Yet constrained by the chemistry
of life, given access to a toolkit of just 20 amino acids, often the same catalytic
patters emerge. This means that by focusing on catalytic sites, functional convergence
around similar mechanisms and structures can be detected. This makes **EnzyMM**
fundamentally knowledge based and interpretable.

While proteins are commonly explored and compared using sequence based methods,
these suffer various drawbacks. Principally it is not immediately apparent which part
of the sequence is most important. Secondly, sequence diverges much faster than
structure obscuring similarities as proteins mutate and evolve. Sequence rearrangements
can completely break certain approaches.

With **PyJess** and **EnzyMM** we develop an extremely fast yet precise way to identify
known catalytic patters without any of the drawbacks of sequence based methods. 
Thanks to very detail-oriented yet large scale curation efforts in the M-CSA and some
`clever engineering <https://pyjess.readthedocs.io/en/latest/guide/optimizations.html>`_
we can scan most protein structures against nearly 7000 catalytic arrangements in tenths
of a second. Perfectly suited for exploring the expanding universe of predicted protein
structures, we are excited to present **EnzyMM**!


🔧 Setup
-----

**EnzyMM** is available for all modern Python versions (3.7+).

Run ``pip install enzymm`` in a shell to download the latest release from PyPI,
or have a look at the :doc:`Installation page <guide/install>` to find other ways 
to install **EnzyMM**.

Alternatively, we provide both Docker and Apptainer Images 🖼️.

📚 Library
-------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   User Guide <guide/index>
   API Reference <api/index>

🔎 Related Projects
----------------

**EnzyMM** draws on a library of catalytic templates derived by Ioannis Riziotis from
the `Mechanism and Catalytic Site Atlas <https://www.ebi.ac.uk/thornton-srv/m-csa/>`_
which is developed at the `Thornton Group <https://www.ebi.ac.uk/research/thornton/>`_
and maintained by Antonio Ribeiro.

**EnzyMM** utilizes **PyJess**, a `Cython <https://cython.org/>`_
wrapper of `JESS <https://github.com/iriziotis/jess>`_,
for finding matches between templates and query structures.
**PyJess** was developed by Martin Larralde.
Please have a look at `PyJess <https://github.com/althonos/pyjess>`_!

🔖 Citations
---------
**EnzyMM** is academic software but relies on many previous approaches.  
**EnzyMM** itself can not yet be cited but a preprint is in preparation.
We intend to publish during the autumn of 2025.  

We kindly ask you to cite both:

- **PyJess**, for instance as:
   - PyJess, a Python library binding to Jess (Barker *et al.*, 2003).

- **Mechanism and Catalytic Site Atlas** as:
   - Ribeiro AJM et al. (2017), Nucleic Acids Res, 46, D618-D623. Mechanism and Catalytic Site Atlas (M-CSA): a database of enzyme reaction mechanisms and active sites. DOI:10.1093/nar/gkx1012. PMID:29106569.

⚖️ License
-------

**EnzyMM** is provided under the `MIT License <https://choosealicense.com/licenses/mit/>`_.
Both **PyJess** and the original **JESS** code are distributed under the 
`MIT License <https://choosealicense.com/licenses/mit/>`_ as well.

Though conceived at the `EMBL-EBI <https://www.ebi.ac.uk/>`_ in Hinxton,
UK in the `Thornton Group <https://www.ebi.ac.uk/research/thornton/>`_, **EnzyMM** is now
developed by Raymund Hackett and the `Zeller Group <https://zellerlab.org/>`_ at
the `Leiden University Medical Center <https://www.lumc.nl/en/>`_ 
in Leiden in the Netherlands with continuing support from the Thornton Group.
