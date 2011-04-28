========================
Tutorial on ClusterBeans
========================

.. testsetup::
   
   import numpy
   numpy.random.seed(1221)
   
This tutorial shows the basic procedure of detecting, sorting, and exporting the spike data 
with the SpikeSort software. Here we use the built-in component based command-line interface
SpikeBeans, which is abstract and very easy to use.
 
You will also find the information on how to adjust the provided scripts for your needs.
 
.. note::
	The procedure described in this tutorial, can also be implemented using the SpikeSort's methods directly (See "Tutorial on manual sorting").

The procedure consists of several steps:

* Reading the data
* Spike sorting
* Exporting results


1. Reading the data
-------------------

Before you can start spike sorting you have to load data with raw extracellular
recordings. Such data can be obtained from microelectrodes, tetrodes or shank
electrodes. SpikeSort currently supports data saved in Bakerlab and HDF5 format
but new formats can be easily added.
   
To start with, you can download a sample data file from 
:download:`here <../../data/tutorial.h5>`.
  
This hierarchical file is organized as following::
   
   /subject/session/electrodeid
      
For example, data from electrode `el1` recorded in session `session01` of 
subject `SubjectA` can be found under `/SubjectA/session01/el1`
   
I will assume that you downloaded this file and saved it to :file:`data` 
directory.

The fastest way prepare the data, is to use the provided script :file:`cluster_beans.py`
which can be found in the :file:`SpikeSort/examples/sorting/` directory. By default, the script
is configured and ready to use with the tutorial data. For now, we'll use the unmodified version.
However, you can find the detailed information on how to adapt the script for your needs in the end
of this tutorial.

Before you run the script, create the environment variable called DATAPATH which points to the folder
where the data is located. On UNIX-based systems, it can be added from the terminal as follows::

	$ export DATAPATH='where/your/data/is/located'
	
Then, just fire it up from the python shell in the interactive mode::

	$ python -i cluster_beans.py

or, alternatively, within the ipython interactive shell::

	%run cluster_beans.py

Now everything is ready to start sorting the spikes.

2. Spike sorting
----------------

1. **Visualizing the data**

   To provide some useful data visualization, the :file:`cluster_beans.py` creates four plotting objects:
	
   * *browser* - the Spike Browser
   * *plot1* - feature space viewer
   * *plot2* - spike waveshape viewer
   * *legend* - no comments
	
   In order to see the output of such an object, just call its :py:func:`show` method, e.g.
   
   .. doctest::
   
   		>>> browser.show()	
   		
   |
   
   *Spike Browser*
   
   `TODO: add the figure here`
   
   The four horizontal black curves are the [filtered] voltage traces recorded
   from different channels (sorted by id from bottom to up) of the electrode
   `el1` (can be changed in the script). The colored pieces are the detected spikes'
   waveshapes. The cell-color correspondence can be found in the legend.
   
   Use the "+" and "-" keys to scale the vertical axis, and the "Prev" and "Next"
   buttons to navigate across the temporal axis.
   
   |
   
   *Feature space viewer*
   
   `TODO: add the figure here`
   
   To sort the spikes, some characteristic features that may be used to differentiate
   between the waveshapes have been calculated (e.g. peak-to-peak amplitude,
   projections on the principal components).
   
   To help the user identify the features, all features are assigned with abbreviated
   labels. For example, feature ``Ch0:P2P`` denotes peak-to-peak amplitude in contact 
   (channel) 0.
   
   The Feature space viewer plots the two-dimensional projections of the feature space
   and histograms of features.
   
   .. note::
   
       Depending on how many features are viewed, the subplots may be too small.
       To zoom in/out the subplot, target it with the mouse and press the "z" key.
       
   |
       
   *Spike waveshape viewer*
   	
   `TODO: add the figure here`
   	
   This component plots the aligned and overlapped spike waveshapes. The spikes
   recorded from different channels are shown in different subplots, numbered in the
   left-right, top-down way. 
   
   You can also zoom the subplots here as in the Feature space viewer.
   
   |
   
   *Legend*
   	
   `TODO: add the figure here`
   
   For the convenience, the legend is plotted on the separate figure with this
   component.
   
   |
   
#. **Managing the spikes**
   
   The aim of the spike sorting is to differentiate one or several cells' firing
   from other unnecessary activity (such as background noise or stimulus artifacts).
   This can be partially done by the automatic clustering in the feature space.
   However, for the reliable results, some manual manipulations are needed and the
   best settings have to be identified using trial-and-error procedure. It usually
   involves removing/merging cells (clusters), reclustering the data, and changing
   the spike detection threshold.
   
   Before we proceed, it will be convenient to create some references:
   
   .. doctest::
   
   	  >>> ca = base.features['LabelSource']         # points to the ClusterAnalyzer instance
   	  >>> sd = base.features['SpikeMarkerSource']   # points to the SpikeDetector instance

   Looking at the spike waveshapes, one might find, that the Cell 3 is most probably
   not really a cell, but some noise. Thus, it is not interesting and we want to
   **remove** it.
  
   To remove one or more cells (i.e. clusters), you have to look up their id's
   in the legend and then pass them as arguments to the :py:func:`ca.delete_cells` function:
   
   .. doctest::
   
      >>> ca.delete_cells(3)
      
   After we got rid of the unnesessaey stuff, the waveshape plot looks as follows:
   
   `TODO: add plot here`
   
   All the deleted cells are now assigned the id 0, which can be considered as a trash.
   
   The clusters 1 and 2 are probably contain responses of the same cell (their mean
   waveshapes are very similar). It will be then very convenient to **merge** them into a
   single cluster.
     
   The merging procedure is similar to deletion:
   
   .. doctest::
   
      >>> ca.merge_cells(1,2)
      
   `TODO: add plot here`
   
   If we now look at the feature space viewer, we will see that the Cell4's spikes seem to form
   two clusters instead of just one. So we better **recluster** it into two new cells.
   
   `DISCUSS WITH BARTOSZ: I don't think the above is applicable here...`

   To do so, use the :py:func:`ca.recluster` function:
   
   .. doctest::
   
      >>> ca.recluster(4, 'gmm', 2)
      
   where the arguments are: `cell to recluster`, `clustering algorithm` [#f1]_, and the necessary
   `number of new clusters`
   
   `TODO: add plot here`

   In practice, it may happen that the **threshold** used during the spike detection is too
   high to detect some important activity or too low to leave the noise out. In this case
   you can easily change it (as well as any other SpikeDetector property) adjusting the
   corresponding SpikeDetector property:
   
   .. doctest::
   
      >>> sd.threshold = 90
      >>> sd.update()
      
   `TODO: MAYBE add plot here`
      
3. Exporting the results
------------------------   
   
   Once you done with the cells' differentiation, it is necessary to save the results
   somewhere. Depending on the type of the data used, the differentiated spike times
   can be stored differently. The tutorial data is in the *HDF5* format, so the
   results will be stored inside the initial :file:`tutorial.h5` file.
   
   To export the data we'll use an instance of the :py:class:`ExportCells` component
   which is already created by :file:`cluster_beans.py`. So make sure the data file
   is writable for python and run:
   
   .. doctest::
   
      >>> export.export()
      
      
   Good luck!!!

   `CLARIFY: this didn't work in my case though...`

Appendix: Configuring the spike_beans.py script
-----------------------------------------------

The example script :file:`spike_beans.py` can be easily adjusted to fit your
needs and used with the real data. Here we list the number of fields you might
want to adjust:

* **hd5file**		is the name of the data file (e.g. `\'tutorial.h5\'`)
* **dataset** 		specifies the data we are interested in (e.g. `/SubjectA/session01/el1`)
* **contact** 		sets the contact (channel) for the initial spike detection (e.g. `3`)
* **type** 			the type of spike waveshapes' alignment (e.g. `\'max\'` - align by the peak value)
* **thresh** 		sets the threshold for the automatic spike detection in milivolts ??? (e.g. `70`)
* **filter_freq** 	specifies the filter properties in the form `(pass freq, cut-off freq)` (e.g. `(800., 100.)`)
* **sp_win** 		specifies the window for spike alignment (e.g. `[-0.6, 0.8]`)

Additionally, you can add some features to be taken into account during clustering
and sorting, using the :py:func:`add_feature` function of the
:py:class:`FeatureExtractor` instance. Again, it's pretty intuitive.

Adding the Peak-to-Peak feature:

.. doctest::

   >>> base.features["FeatureSource"].add_feature("P2P")
   
Adding 3 Principal Components to the feature list:

.. doctest::

   >>> base.features["FeatureSource"].add_feature("PCs", ncomps=3)

|
|
|
   
.. [#f1] There are several automatic, semi-automatic and manual methods for clustering.
   Their performance and accuracy depends to large degree on a particular dataset
   and recording setup. In SpikeSort you can choose from several available methods,
   whose names are given as the first argument of :py:func:`spike_sort.cluster.cluster`
   method. The 'gmm' shortcut used in this example, means the Gaussian Mixture Model algorithm
