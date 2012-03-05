.. testsetup::
   
   import numpy
   numpy.random.seed(1221)


Data Structures
===============

.. _raw_recording:

To achieve best compatibility with external libraries most of the data
structures are standard Python *dictionaries*, with at least one key -- `data`. The
`data` key contains the actual data in an array-like object (such as
NumPy array). Other attributes provide metadata that are required by
some methods.


Raw recording
-------------

Raw electrophysiological data sampled at equally spaced time points. It
can contain multiple channels, but all of them need to be of the same
sampling frequency and duration (for example, multiple contacts of
a tetrode). 

The following keys are defined:

  :data: *array*, required 
      
      array-like object (for example :py:class:`numpy.ndarray`) of 
      dimensions (N_channels, N_samples)
  
  :FS: *int*, required
      
      sampling frequency in Hz
  
  :n_contacts: *int*, required 
      
      number of channels (tetrode contacts). It is equal to the size
      of the first dimension of `data`.

.. note::
   
   You may read/write the data with your own functions, but to make the
   interface with the SpikeSort a bit cleaner, you might also want to
   define your custom IO filters (see :ref:`io_filters`)

.. rubric:: Example

We will read the raw tetrode data from :ref:`tutorial_data` using
standard :py:class:`~spike_sort.io.filters.PyTablesFilter`:

  >>> from spike_sort.io.filters import PyTablesFilter
  >>> io_filter = PyTablesFilter('../data/tutorial.h5')
  >>> raw_data = io_filter.read_sp('/SubjectA/session01/el1')
  >>> print(raw_data.keys())          # print all keys
  ['n_contacts', 'FS', 'data']
  >>> shape = raw_data['data'].shape # check size
  >>> print "{0} channels, {1} samples".format(*shape)    
  4 channels, 23512500 samples
  >>> print(raw_data['FS'])    # check sampling frequency
  25000


.. _spike_times:

Spike times
-----------

A sequence of (sorted) time readings at which spikes were generated
(or other discrete events happened). 

The data is store in a dictionary with following keys:

  :data: *array*, required
      
      one-dimensional array-like object with event times (in
      milliseconds) 

  :is_valid: *array*, optional

      boolean area of the same size as `data` -- if an element is False
      the event of the same index is masked (or invalid)

.. note::
   
   You may read/write the data with your own functions, but to make the
   interface with the SpikeSort a bit cleaner, you might also want to
   define your custom IO filters (see :ref:`io_filters`)

.. rubric:: Example

:py:func:`spike_sort.core.extract.detect_spikes` is one of functions
which takes the raw recordings and returns spike times dictionary:


   >>> import numpy as np
   >>> raw_dict = {
   ...       'data': np.array([[0,1,0,0,0,1]]),
   ...       'FS'  : 10000,
   ...       'n_contacts': 1
   ...            }
   >>> from spike_sort.core.extract import detect_spikes
   >>> spt_dict = detect_spikes(raw_dict, thresh=0.8)
   >>> print(spt_dict.keys())
   ['thresh', 'contact', 'data']
   >>> print('Spike times (ms): {0}'.format(spt_dict['data']))
   Spike times (ms): [ 0.   0.4]


Note that in addition to the required data key, 
:py:func:`~spike_sort.core.extract.detect_spikes`
appends some extrcontact a attributes: :py:attr:`thresh` (detection threshold)
and :py:attr:`contact` (contact on which spikes were detected). These
attributes are ignored by other methods.

.. _spike_wave:

Spike waveforms
---------------

Spike waveform structure contains waveforms of extracted spikes. It may be
any mapping data structure (usually a dictionary) with following keys:

:data: *array*, required

    three-dimensional array-like object of size (N_points, N_spikes,
    N_contacts), where:
     
      * `N_points` --  the number of data points in a single waveform, 
      * `N_spikes` -- the total number of spikes and
      * `N_contacts` -- the number of independent channels (for example 4 in a
        tetrode) 

:time: *array*, required

    Timeline of the spike waveshapes (in miliseconds). It must be of
    the same size as the first dimension of data (`n_pts`).

:FS: *int*, optional
     
    Sampling frequency.


:n_contacts: *int*, optional

   Number of independent channels with spike
   waveshapes (see also :ref:`raw_recording`).
  
:is_valid: *array*, optional

   boolean area of the size of second dimension of `data` (N_spikes) -- if an element is False
   the spike with the same index is masked (or invalid)


.. rubric:: Example

Spike waveforms can be extracted from raw recordings (see :ref:`raw_recording`)
given a sequence of spike times (see :ref:`spike_times`) by means of
:py:func:`spike_sort.core.extract.extract_spikes` function:

   >>> from spike_sort.core.extract import extract_spikes
   >>> raw_dict = {     
   ...             'data': np.array([[0,1,1,0,0,0,1,-1,0,0, 0]]),
   ...             'FS'  : 10000,
   ...             'n_contacts': 1
   ...            }      # raw signal
   >>> spt_dict = {      
   ...             'data': np.array([0.15, 0.65, 1])}
   ...            }      # timestamps of three spikes
   >>> sp_win = [0, 0.4] # window in which spikes should be extracted
   >>> waves_dict = extract_spikes(raw_dict, spt_dict, sp_win)

Now let us investigate the returned spike waveforms structure:

* keys:

   >>> print waves_dict.keys()
   ['is_valid', 'FS', 'data', 'time']

* data array shape:

   >>> print(waves_dict['data'].shape)    
   (4, 3, 1)

* extracted spikes:

   >>> print(waves_dict['data'][:,:,0].T) # data contains three spikes
   [[ 1.  1.  0.  0.]
    [ 1. -1.  0.  0.]
    [ 0.  0.  0.  0.]]
   >>> print(waves_dict['time'])          # defined over 4 time points
   [ 0.   0.1  0.2  0.3]

* and potential invalid (truncated spikes):

   >>> print(waves_dict['is_valid'])     # last spike is invalid (truncated)
   [ True True False]

Note that the :py:attr:`is_valid` element of truncated spike is
:py:data:`False`.

.. _spike_features:

Spike features
--------------

This data structure contains features calculated from spike waveforms
using one of the methods defined in :py:mod:`spike_sort.core.features` module 
(one of the :py:func:`fet*` functions, see :ref:`features_doc`). 

The spike features dictionary consits of following keys:

:data: *array*, required
      
    two-dimensional array of size (N_spikes, N_features) that contains
    the actual feature values

:names: *list of str*, required

    list of length N_features containing feature labels  

:is_valid: *array*, optional

   boolean area of of length N_spikes; if an element is False
   the spike with the same index is masked (or invalid, see also
   :ref:`spike_wave`)


.. rubric:: Example

Let us try to calculate peak-to-peak amplitude from some spikes
extracted in :ref:`spike_wave`:

   >>> from spike_sort.core.features import fetP2P
   >>> print(waves_dict['data'].shape) # 3 spikes, 4 data points each
   (4, 3, 1)
   >>> feature_dict = fetP2P(waves_dict)
   >>> print(feature_dict.keys())
   ['is_valid', 'data', 'names']
   >>> print(feature_dict['data'].shape)
   (3, 1)

Then we have one feature for 3 spikes. Let check whether the peak-to-peak amplitudes 
are correctly calculated:
   
   >>> print(feature_dict['data'])
   [[ 1.]
    [ 2.]
    [ 0.]]

as expected (compare with example above). There is only one
peak-to-peak (`P2P`) feature on a single channel (`Ch0`) and its name
is:
   
   >>> print(feature_dict['names'])
   ['Ch0:P2P']

The mask array is inherited from :py:data:`waves_dict`:

   >>> print(feature_dict['is_valid'])
   [ True True False]

.. _spike_labels:

Spike labels
------------

Spike labels are the identifiers of a cell (unit) each spike was
classified to. Spike labels are **not** dictionaries, but arrays of
integers -- one cluster index per spike.

.. rubric:: Example

Let us try to cluster the spikes described by `Sample` feature using
K-means with K=2:
   
   >>> from spike_sort.core.cluster import cluster
   >>> feature_dict = {
   ...                  'data' : np.array([[1],[-1], [1]]),
   ...                  'names' : ['Sample']
   ...                }
   >>> labels = cluster('k_means', feature_dict, 2)
   >>> print(labels)
   [1 0 1]

As expected :py:data:`labels` is an array describing two clusters: 0 and 1.
        
