Cluster (:mod:`spike_sort.cluster`)
===================================

.. currentmodule:: spike_sort.core.cluster

Module with clustering algorithms.

Utility functions
-----------------

Spike sorting is usually done with the
:py:func:`~spike_sort.core.cluster.cluster` function which takes as an
argument one of the clustering methods (as a string).

Others functions help to manipulate the results:

.. autosummary:: 
   
   cluster
   split_cells



Clustering methods
------------------

Several different clustering methods are defined in the module. Each
method should take at least one argument -- the features structure.

.. autosummary::

   k_means_plus
   gmm
   manual
   none
   k_means
   

Reference
---------

.. automodule:: spike_sort.core.cluster
   :members:
