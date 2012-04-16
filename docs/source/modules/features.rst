Calculate features (:mod:`spike_sort.features`)
===============================================

.. currentmodule:: spike_sort.core.features 

Provides functions to calculate spike waveforms features.

.. _features_doc:

Features
--------

Functions starting with `fet` implement various features calculated
from the spike waveshapes. They have usually one required argument
:ref:`spike_wave` structure and some can have optional arguments (see
below)

Each of the function returns a mapping structure (dictionary) with the following keys:

 * `data` -- an array of shape (n_spikes x n_features)
 * `names` -- a list of length n_features with feature labels

The following features are implemented:

   .. autosummary:: 
   
      fetPCs
      fetP2P
      fetSpIdx
      fetSpTime
      fetSpProjection

Tools
-----

This module provides a few tools to facilitate working with features
data structure:

.. autosummary::
   
   split_cells
   select
   combine
   normalize

Auxiliary
---------
.. autosummary:: 

   PCA
   add_mask


Reference
---------

.. automodule:: spike_sort.core.features
   :members:
