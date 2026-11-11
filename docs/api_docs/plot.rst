.. _api_plot:

****
Plot
****

Counts
======

.. autoclass:: moonstone.plot.counts.PlotCountsStats
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:

.. autoclass:: moonstone.plot.counts.PlotTaxonomyCounts
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:

Metadata
========

.. autoclass:: moonstone.plot.metadata.PlotMetadataStats
    :members:
    :undoc-members:
    :special-members: __init__
    :inherited-members:
    :member-order: bysource

Category flows
==============

See :ref:`category_flows`.

.. autofunction:: moonstone.plot.category_flows.plot_category_flows

.. autofunction:: moonstone.plot.category_flows.plot_taxonomy_flows

.. autoclass:: moonstone.plot.graphs.sankey.SankeyDiagram
    :members: plot_one_graph
    :special-members: __init__

.. autoclass:: moonstone.plot.graphs.parcats.ParallelCategoriesGraph
    :members: plot_one_graph
    :special-members: __init__
