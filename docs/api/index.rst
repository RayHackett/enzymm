API Reference
=============

.. currentmodule:: enzymm

.. automodule:: enzymm

Input
-----

.. autosummary::
    :nosignatures:

    jess_run.load_molecules
    template.load_templates

.. toctree::
    :maxdepth: 1
    :hidden:

    input

Templates
---------

.. autosummary::
    :nosignatures:

    template.load_templates
    template.Vec3
    template.Residue
    template.AnnotatedResidue
    template.Cluster
    template.Template
    template.AnnotatedTemplate

.. toctree::
    :maxdepth: 1
    :hidden:

    templates


Matching
--------

.. autosummary::
    :nosignatures:

    jess_run.Matcher
    jess_run.Match
    jess_run.ModelEnsemble
    jess_run.LogisticRegressionModel

.. toctree::
    :maxdepth: 1
    :hidden:

    matching

M-CSA
-----

.. autosummary::
    :nosignatures:

    mcsa_info.HomologousPDB
    mcsa_info.HomologousResidue
    mcsa_info.NonReferenceCatalyticResidue
    mcsa_info.ReferenceCatalyticResidue

.. toctree::
    :maxdepth: 1
    :hidden:

    m-csa

Output
------

.. autosummary::
    :nosignatures:

    output.BaseTable
    output.FullMatchTable
    output.SimpleMatchTable
    output.MatchResidueTable
    output.Tables

.. toctree::
    :maxdepth: 1
    :hidden:

    output