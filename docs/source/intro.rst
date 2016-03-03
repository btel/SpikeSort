Introduction
============

What is spike sorting?
----------------------

    *Spike sorting is a class of techniques used in the analysis of
    electrophysiological data. Spike sorting algorithms use the
    shape(s) of waveforms collected with one or more electrodes in the
    brain to distinguish the activity of one or more neurons from
    background electrical noise.*

         Source: Spike sorting, Wikipedia_

Spike sorting usually consists of the following steps [Quiroga2004]_:

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


#. Sorting evaluation (*evaluate*)

   After sorting is done and we have determined spike times of
   different units, we have to make sure that the quality of sorting
   is sufficient for further analysis. There are multiple visual and
   statistical methods that help in this evaluations, such as
   inter-spike intervals histograms, waveforms overlays,
   cross-corellograms etc. [Hill2011]_.

SpikeSort is a comprehensive library whose goal is to accompany you trough
the entire process of spike sorting - from detection to evaluation.
It is not fully automatic or on-line sorting program. Although we
include lots of functions, that help to automatize some of the
repetitive task, good sorting will always require human supervision.

Design goals
------------

There are several of spike sorting programs on the market both commercial and free
and open source. We started this projects to offer an alternative that
would be free, scriptable and preferably in our favourite language:
Python. However, we did not try to re-invent the
wheel and whenever we could we leveraged the established libraries.

We had a few design goals when we worked on SpikeSort;

* modular

  Spike sorting is modular - there are several steps and different
  techniques that can and
  should be mixed-and-matched to adjust the process to the data we
  try to analyze. Spike sorting library should be modular as well to allow for
  much flexibility. This is achieved by composing the library of
  independent components that can be easily inserted into or deleted
  from the spike sorting workflow.

* customizable

  No two data sets are the same. Different experimental protocols,
  different acquisition systems, different neural systems result in
  different properties of the data and thus require different methods.
  Therefore, a good spike sorting library should allow for easy and
  flexible customizations of the algorithms used at each stage of spike
  sorting and their parameters.

* easy-to-use

  Flexibility usually comes at price: the usability. Nevertheless, a
  transparent design can allow complex systems to be user-friendly.
  The interface should allow to focus on the data and not on the
  peculiarities of specific software solutions. 

* fast

  In practice, one will try to discriminate thousands of spikes of
  tens of different neurons. Any performance optimizations will save
  you precious minutes (hours or even days) for the task you are most
  interested in: decoding what the cells actually do.

* compatible with standard libraries (NumPy, SciPy, matplotlib)

  There is no need to reinvent the wheel. Python developers provide
  hundreds of optimized, well-tested and widely-used libraries with
  great community support. Why not use them? Moreover, any data coming
  from the spike sorting libraries should be easily pluggable into
  third-party analysis routines and databases.

Although still much work is required to meet all the goals, we kept
them all in mind while designing the SpikeSort. As a result, SpikeSort
is already  an usable and powerful framework that will help you to
get most of your data.

Why Python?
-----------

Python is a interpreted and very dynamic language with huge very
enthusiastic community. It grew to be the de-facto standard in 
computational science [Langtangen2009]_ and it rapidly gains momentum in
experimental disciplines [Hanke2009]_. Python is also completely free and available
for multiple platforms - it means that you can run it at home, at work
or give it to your students with no additional costs. Last but not
least Python can be easily interfaced with other languages making it a
''glueing'' language that can make two independent libraries to
communicate. 

Installation
------------

You can download the most recent release of SpikeSort from github::

   git clone git://github.com/btel/SpikeSort.git

In order to install SortSpike you need following libraries:

* python 2.6 or 2.7
* setuptools
* scipy
* numpy
* pytables
* matplotlib (only for plotting)

Optional dependencies are:

* scikits.learn - clustering algorithms
* neurotools - spike train analysis
* ipython - enhanced python shell

If some of the python packages are not available on your system you
can install them with easy-install::

   easy_install numpy scipy pytables matplotlib

.. note::

   If you are not familiar with Python packaging system we recommend
   you installing a complete Python distribution from a company called
   Enthought: `EPD <http://www.enthought.com/products/epd.php>`_
   (there are free academic licenses). Installers for Windows, MacOSX
   and Linux are available.

If you have the above libraries you can install SpikeSort simply
issuing the command::

   python setup.py install

If you prefer to install it in your home directory you may try::

   python setup.py install --user

but remember to add :file:`$HOME/.local/lib/python2.6/site-packages` to your python
path.

After a successful installation you can run the supplied tests::

   python setup.py nosetests

If you don't have all the optional dependencies, be prepared for some
tests errors.

Examples
--------

In :file:`examples/sorting` subdirectory you will find some sample scripts,
which use SpikeSort for spike sorting

* :file:`cluster_manual.py` - sort spikes by manual cluster cutting
* :file:`cluster_auto.py` - automatically cluster with GMM (Gaussian
  Mixture Models) algorithm (see our tutorial
  :ref:`lowlevel_tutorial`)
* :file:`cluster_beans.py` - run full stack spike-sorting evnvironment
  and show spikes in a spike browser (see our tutorial :ref:`beans_tutorial`)

In order to run these examples, you need to download :ref:`tutorial_data` and define an environment variable ``DATAPATH``::

   export DATAPATH=/path/to/data/directory

where ``/path/to/data/directory`` points to the directory where you
downloaded the data file.

Once you have the tutorial data, you may run above script, for example::

   python -i cluster_auto.py

.. note::

   The ``-i`` in python command will leave a Python interpreter open
   for interactive exploration - read more  in our tutorials.

Similar software
----------------

There are a few open source packages for spike sorting that have
different design and use case:

* `spyke <http://www.frontiersin.org/Neuroinformatics/10.3389/neuro.11.009.2008/pdf/abstract>`_

* `OpenElectrophy <http://packages.python.org/OpenElectrophy/>`_

* `spikepy <http://code.google.com/p/spikepy/>`_ 
  
* `Klusters <http://klusters.sourceforge.net/>`_ 




References
----------

.. _tutorial.h5: https://github.com/btel/SpikeSort/releases/download/v0.12/tutorial.h5

.. _Wikipedia: http://en.wikipedia.org/wiki/Spike_sorting

.. [Quiroga2004] Quiroga, RQ, Z. Nadasdy, Y. Ben-Shaul, and others. *“Unsupervised Spike Detection and Sorting with Wavelets and Superparamagnetic Clustering.”* Neural Computation **16**, no. 8 (2004):1661. `<http://www.vis.caltech.edu/~rodri/papers/Spike_sorting.pdf>`_

.. [Hill2011] Hill, Daniel N, Samar B Mehta, and David Kleinfeld. *“Quality Metrics to Accompany Spike Sorting of Extracellular Signals.”* The Journal of Neuroscience **31**, no. 24 (2011): 8699-8705. `<http://neuro.cjb.net/content/31/24/8699.full>`_

.. [Hanke2009] Hanke, Michael, Yaroslav O. Halchenko, Per B. Sederberg, Emanuele Olivetti, Ingo Fründ, Jochem W. Rieger, Christoph S. Herrmann, James V. Haxby, Stephen José Hanson, and Stefan Pollmann. *“PyMVPA: a Unifying Approach to the Analysis of  Neuroscientific Data.”* Frontiers in Neuroinformatics **3** (2009): 3.

.. [Langtangen2009] Langtangen, Hans Petter. *Python Scripting for Computational Science*. 3rd ed. Springer, 2009.

