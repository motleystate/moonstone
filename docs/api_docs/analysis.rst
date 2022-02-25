.. _api_analysis:

********
Analysis
********

Differential Analysis
=====================

.. currentmodule:: moonstone.analysis.differential_analysis

.. autosummary::
    :toctree: stubs

    DifferentialAnalysis


Diversity
=========

Alpha-Diversity
"""""""""""""""

Alpha-diversity corresponds to the intra-sample diversity

.. autoclass:: moonstone.analysis.diversity.alpha.ShannonIndex
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:

.. autosummary::
    :toctree: stubs

    moonstone.analysis.diversity.alpha.SimpsonInverseIndex
    moonstone.analysis.diversity.alpha.Chao1Index
    
Phylogenetic alpha-diversity

.. autoclass:: moonstone.analysis.diversity.alpha.FaithsPhylogeneticDiversity
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:

==================================================================================================

Beta-Diversity
"""""""""""""""

Beta-diversity corresponds to the inter-samples diversity

.. autoclass:: moonstone.analysis.diversity.beta.BrayCurtis
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:

Phylogenetic beta-diversity

.. autoclass:: moonstone.analysis.diversity.beta.WeightedUniFrac
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:

.. autosummary::
    :toctree: stubs

    moonstone.analysis.diversity.beta.UnweightedUniFrac
