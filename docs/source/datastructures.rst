Data Structures
===============

.. _raw_recording:

Raw recording
-------------

.. _spike_times:

Spike times
-----------

.. _spike_wave:

Spike waveforms
---------------

Spikewave structure contains waveforms of extracted spikes. It may be
any mapping datastructure (usually a dictionary) with following keys:

.. py:attribute:: data
  
  3D array (numpy or pytables) of shape `(n_pts,
  n_spikes, n_contacts)`, where 
  
  * `n_pts` -- number of data points in a single waveform, 
  
  * `n_spikes` -- the total number of spikes and
  
  * `n_contact` -- the number of independent channels (for example 4 in a
    tetrode) 

.. py:attribute:: time

   Timeline of the spike waveshapes (in miliseconds). It must be of
   the same size as the first dimension of data (`n_pts`).

.. py:attribute:: FS

   Sampling frequency.


.. py:attribute:: n_contacts

   (Auxilarly variable) Number of independent channels with spike
   waveshapes (see also :py:attr::`data` definition).

.. _spike_features:

Spike features
--------------

.. _spike_labels:

Spike labels
------------


