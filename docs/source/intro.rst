Introduction
============

What is spike sorting?
----------------------

1. Spike detection (*detect*)
#. Spike waveform extraction (*extract*)
#. Feature extraction (*feature*)
#. Clustering (*cluster*)
#. Sorting evaluation (*evalulate*)

Desgin goals
------------

* extensible
* easy-to-use
* fast
* compatible with other libraries

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

