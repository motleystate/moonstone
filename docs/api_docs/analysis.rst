.. _api_analysis:

********
Analysis
********

Differential Analysis
=====================

.. automodule:: moonstone.analysis.differential_analysis
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:


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

.. autoclass:: moonstone.analysis.diversity.alpha.SimpsonInverseIndex
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:

.. autoclass:: moonstone.analysis.diversity.alpha.Chao1Index
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:
    
Phylogenetic alpha-diversity includes the phylogenetic closeness of species in the computation of the diversity indexes

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

Phylogenetic beta-diversity includes the phylogenetic closeness of species in the computation of the diversity indexes

.. autoclass:: moonstone.analysis.diversity.beta.WeightedUniFrac
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:

.. autosummary::
    :toctree: stubs

    moonstone.analysis.diversity.beta.UnweightedUniFrac
