Introduction
============

What is spike sorting?
----------------------

    *Spike sorting is a class of techniques used in the analysis of
    electrophysiological data. Spike sorting algorithms use the
    shape(s) of waveforms collected with one or more electrodes in the
    brain to distinguish the activity of one or more neurons from
    background electrical noise.*

         Source, Wikipedia_

Spike sorting usually consists of the following steps:

1. Spike detection (*detect*)

   Spikes are very rapid and often sparse events, so that they appear
   only a tiny fraction of the recordings. To achieve a sort of
   compression and easy the subsequent analysis, the times of
   occurrence of putative spikes are first identified in the continuous
   by means of thresholding. This method return only events that cross
   a specified threshold (selected based on some data statistics, such
   as standard deviation, or visually).

#. Spike waveform extraction (*extract*)

   Based on the series of spike times identified in the previous step,
   we now may proceed to extract the spike waveforms by taking a small
   segment of the signal around each spike time. Such segments are
   usually automatically aligned to a specific feature of the
   waveform, for example maximum or minimum.

#. Feature extraction (*feature*)

   This is one of the most important steps in which the silent
   features of the spikes such as peak-to-peak amplitude or spike
   width are calculated based on spike waveshapes. The features should
   be preferably low-dimensional and should well differentiate spikes
   of different cells and noise.

#. Clustering (*cluster*)

   At the heart of spike sorting is the clustering that uses
   automatic, semi-automatic or manual methods to identify groups of
   spikes belonging to the same cell (often called unit). The
   procedure is usually applied in n-dimensional space of spike
   features, where each feature is a single dimension. Since it is
   very difficult to do that visually for more than 2 features (on a
   plane), one usually resorts to different clustring algorithms, such
   as K-means or Gaussian Mixture Models that can handle large number
   of dimension.


#. Sorting evaluation (*evalulate*)

   After sorting is done and we have determined spike times of
   different units, we have to make sure that the quality of sorting
   is sufficient for further analysis. There are multiple visual and
   statistical methods that help in this evaluations, such as
   inter-spike intervals histograms, waveforms overlays,
   cross-corellograms etc.

SpikeSort is a comprehensive library whose goal is to accompany you trough
the entire process of spike sorting - from detection to evaluation.
It is not fully automatic or on-line sorting program. Although we
include lots of functions, that help to automatize some of the
repetitive task, good sorting will always require human supervision.

.. _Wikipedia: http://en.wikipedia.org/wiki/Spike_sorting


Design goals
------------

There are several of spike sorting programs on the market both commercial and free
and open source. We started this projects to offer an alternative that
would be free and scriptable, perferrably in our favourite language
Python. However, we did not try to re-invent the
wheel and whenever we could we leveraged the established Python libraries.

We had a few design goals when we worked on SpikSort;

* modular

  Spike sorting is modular - there are several steps and different
  techniques that can and
  should be mixed-and-matched to adjust the process to the data we
  try to analyze. Spike sorting library should be modular as well to allow for
  much flexibility.

* customizable

  No two data sets are the same. Different experimental protocols,
  different acquistion systems, different neural systems require
  different methods. Therefore, a good spike sorting library should
  allow for easy and flexible customizations.

* easy-to-use

  Flexibility usually comes at price: the usability. Nevertheless, a
  transparent design can allow complex systems to be user-friendly.
  The interface should allow to focus on the data and not on the
  peculiarities of specific software solutions (we should not be
  limited by the medium).

* fast

  In practice, one will try to discriminate thousends of spikes of
  tens of different neurons. Any performance optimizations will save
  you precious minutes (hours or even days) for the task you are most
  interested in: decoding what the cells actually do.

* compatible with standard libraries (NumPy, SciPy, matplotlib)

  There is no need to reinvent the wheel. Python developers provide
  hundreds of optimized, well-tested and widely-used libraries with
  great community support. Why not use them? Moreover, any data coming
  from the spike sorting libraries should be easily pluggable to
  third-party analysis routines and databases.


Why Python?
-----------

Python is a interpreted and very dynamic language with huge very
enthusiastic community. It grew to be the de-facto standard in 
computational science and it rapidly gains momentum in
experimental disciplines. Python is also completely free and available
for multiple platforms - it means that you can run it at home, at work
or give it to your students with no additional costs. Last but not
least Python can be easily interfaced with other languages making it a
''glueing'' language that can make two independent libraries to
communicate. 

Installation
------------

In order to install SortSpike you need following libraries:

* python 2.6
* setuptools
* scipy
* numpy
* pytables
* matplotlib (only for plotting)

If you have the above libraries you can install SpikeSort simply
issuing the command::

   python setup.py install

If you prefer to install it in your home directory you may try::

   python setup.py install --user

but remember to add :file:`$HOME/.local/lib/python2.6/site-packages` to your python
path.

After a successful installation you can run the supplied tests::

   python setup.py nosetests

Examples
--------

In :file:`examples` subdirectory you will find some sample scripts,
which use SpikeSort:

* :file:`examples/cluster_spikes.py`
* :file:`examples/plot_features.py`

