==========================
Tutorial on manual sorting
==========================

.. testsetup::
   
   import numpy
   numpy.random.seed(1221)

1. **Read data**.

   Before you can start spike sorting you have to load data with raw extracellular
   recordings. Such data can be obtained from microelectrodes, tetrodes or shank
   electrodes. SpikeSort currently supports data saved in Bakerlab and HDF5 format
   but new formats can be easily added.
   
   To start with, you can download a sample data file from 
   :download:`here <../../data/tutorial.h5>`.
   
   This hierachical file is organised as following::
   
      /subject/session/electrodeid
      
   For example, data from electrode `el1` recorded in session `session01` of 
   subject `SubjectA` can be found under `/SubjectA/session01/el1`
   
   I will assume that you downloaded this file and saved it to :file:`data` 
   directory.
   
   You can load this file using one of I/O fiters from 
   :py:mod:`spike_sort.io.filter` module:
   
   .. doctest::
   
      >>> from spike_sort.io.filters import PyTablesFilter
      >>> dataset = '/SubjectA/session01/el1'
      >>> io_filter = PyTablesFilter('../data/tutorial.h5')
      >>> raw = io_filter.read_sp(dataset)
      
   :py:data:`raw` is a dictionary which contains the raw data (in this case it is
   a pytables compressed array) under :py:attr:`data`
   key:

   .. doctest::
   
      >>> print raw['data']
      /SubjectA/session01/el1/raw (CArray(4, 23512500)) ''
      
   The size of the data is 23512500 samples in 4 independent channels (`contacts`
   in the tetrode).
   

#. **Detect spikes**.

   The first step of spike sorting is spike detection. It is usually done by 
   thresholding the raw recordings. Let us use an automatic threshold on 
   contact 3 (remember that indexing always starts with 0!):
   
   .. doctest::
   
      >>> from spike_sort import extract
      >>> spt = extract.detect_spikes(raw,  contact=3, thresh='auto')
      
   Let us see now how many events were detected:
   
   .. doctest::
   
      >>> print len(spt['data'])
      16293
      
   We should make sure that all events are aligned to the same point of reference,
   for example, the maximum amplitude. To this end we first define a window
   around which spikes should be centered and then recalculate aligned event times:
   
   .. doctest::
      
      >>> sp_win = [-0.2, 0.8]
      >>> spt = extract.align_spikes(raw, spt, sp_win, type="max", 
      ...                            resample=10)
      
   `resample` option is optional: it results in upsampling (10 times) the original 
   waveforms to obtain better resolution of event times.
   
   After spike detection and alignment we can finally extract the spike waveforms:
   
   .. doctest::
  
      >>> sp_waves = extract.extract_spikes(raw, spt, sp_win)
      
   The resulting structure is a dictionary whose :py:attr:`data` key is an array
   containing the spike waveshapes. Note that the array is three-dimensional and
   sizes of its dimensions reflect:
     
     * 1st dimmension: number of samples in each waveform,
     * 2nd: number of spikes,
     * 3rd: number of contacts
   
   .. doctest::
   
      >>> print sp_waves['data'].shape
      (25, 15541, 4)
      
   In practice, you do not to take care of such details. However, it is always
   a good idea to take a look at the obtained waveforms. 
   :py:mod:`spike_sort.ui.plotting` module contains various functions which will
   help you to visualize the data. To plot waveshapes you can use: 
   :py:func:``plot_spikes``.
   
   .. doctest::
   
      >>> from spike_sort.ui import plotting
      >>> plotting.plot_spikes(sp_waves, n_spikes=200)
   
   .. plot:: source/pyplots/tutorial_spikes.py
   
   It is apparent from the plot that the spike waveforms of a few different cells
   and also some artifacts were detected. In order to separate these activities, 
   in the next step we will perform *spike sorting*.

#. **Calculate features**.

   Before we can sort spikes, we should calculate some characteristic features 
   that may be used to differentiate between the waveshapes. Module 
   :py:mod:`spike_sort.features` defines several of such features, for example
   peak-to-peak amplitude (:py:func:`fetP2P`) and projections on principal 
   components (:py:func:`fetPCs`). Now, we will calculate peak-to-peak amplitudes
   and PC projections on each of the contact, and then combine them into a single
   object:
   
   .. doctest::
   
      >>> from spike_sort import features
      >>> sp_feats = features.combine(
      ...      (
      ...       features.fetP2P(sp_waves),
      ...       features.fetPCs(sp_waves)
      ...      )
      ... )
   
   To help the user identify the features, all features are assigned with abbreviated
   labels:
   
   .. doctest::
   
      >>> print sp_feats['names']
      ['Ch0:P2P' 'Ch1:P2P' 'Ch2:P2P' 'Ch3:P2P' 'Ch0:PC0' 'Ch1:PC0' 'Ch2:PC0'
       'Ch3:PC0' 'Ch0:PC1' 'Ch1:PC1' 'Ch2:PC1' 'Ch3:PC1']
      
   For examples feature ``Ch0:P2P`` denotes peak-to-peak amplitude in contact 
   (channel) 0.
   
   Let us plot the two-dimensional 
   projections of the feature space and histograms of features:
   
   .. doctest::
  
      >>> plotting.plot_features(sp_feats)
      
   .. plot:: source/pyplots/tutorial_features.py

#. **Cluster spikes**.

   Finally, based on the calculated features we can perform spike clustering. This
   step is a little bit more complex and the best settings have to be identified
   using trial-and-error procedure.
   
   There are several automatic, semi-automatic and manual methods for clustering.
   They performance and accuracy depends to large degree on a particular dataset
   and recording setup. In SpikeSort you can choose from several available methods,
   whose names are given as the first argument of :py:func:`spike_sort.cluster.cluster`
   method.
   
   We will start with an automatic clustering :py:func:`gmm`, which requires
   only the feature object :py:data:`sp_feats` and number of clusters to identify.
   It attempts to find a mixture of gaussian functions which best approximates the
   distribution of spike feature datapoints (gaussian mixture model).
   Since we do not know, how many cells were picked up by the electrode we guess
   an initial number of clusters, which we can modify later on:
   
   .. doctest::
      
      >>> from spike_sort import cluster
      >>> clust_idx = cluster.cluster("gmm",sp_feats,4)
      
   The resulting data is just assigning a number (cluster index) to each spike from
   the feature array :py:data:`sp_feats`.
   
   You can use the plotting module to draw the 
   feature vectors with color reflecting group to which each spike was assigned:
   
   .. doctest::
   
      >>> plotting.plot_features(sp_feats, clust_idx)
      
   .. plot:: source/pyplots/tutorial_clusters.py

   or you can see the spike waveshapes:
   
   .. doctest::
     
      >>> plotting.plot_spikes(sp_waves, clust_idx, n_spikes=200)
      >>> plotting.show()

   .. plot:: source/pyplots/tutorial_cells.py
      
   If you are not satisfied with the results or you think you might do better, 
   you can also try manual sorting using cluster cutting method::
   
      >>> from spike_sort.ui import manual_sort
      >>> cluster_idx = manual_sort.show(features, sp_waves,
      ...                                ['Ch0:P2P','Ch3:P2P'],
      ...                                show_spikes=True)
      
   This function will open a window in which you can draw clusters of arbitrary
   shapes, but beware: you can draw only on two dimensional plane, so that you 
   are limited to only two features!

#. **Export data**.

   Once you are done with spike sorting, you can export the results to a file.
   To this end you can use the same :py:mod:`spike_sort.io` module we used 
   for reading. Here, we will save the spike times of a selected cell
   back to the file we read the data from. 
   
   First, we need to extract the spike times 
   of the discriminated cells:
   
   .. doctest:: 
  
      >>> spt_clust = cluster.split_cells(spt, clust_idx)

   It will create a dictionary whose keys are the cell labels pointing
   to spike times of the specific cell. For example, to extract spike
   times of cell 0:

   .. doctest::

      >>> print spt_clust[0]
      {'data': array([  2.11884000e+02,   2.37192000e+02,   3.45244000e+02, ...,
               9.36228740e+05,   9.36269656e+05,   9.36527580e+05])}
 
      
   Then we may export them to the datafile:

   .. doctest::
   
      >>> from spike_sort.io import export
      >>> cell_template = dataset + '/cell{cell_id}'
      >>> export.export_cells(io_filter, cell_template, spt_clust, overwrite=True)
      
   This will create a new node in :file:`tutorial.h5` containing  spike times of 
   the discriminated cell ``/SubjectA/session01/el1/cell{1-4}``, 
   which you can use for further analysis.
  
   Don not forget to close the I/O filter at the end of your analysis:

   ..doctest::

      >>> io_filter.close()
   
   Good luck!!!
   
   TODO: this must be automated: rewrite io module to provide basic I/O functions
   (read_spt, write_spt, etc.) and implement a module with abstract operations,
   such as write clustering results to a file (''template design pattern'').
   
   
