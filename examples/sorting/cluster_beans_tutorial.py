#############################################
# Adjust these fields for your needs

#data source
hd5file = 'tutorial.h5'
dataset = "/SubjectA/session01/el1"

#spike detection/extraction properties
contact = 3
type = "max"
thresh = 70
filter_freq= (800., 100.)

#conf_file = "../data/gollum.inf"
sp_win = [-0.6, 0.8]

#############################################

import matplotlib
matplotlib.use("TkAgg")
matplotlib.interactive(True)

from spike_beans import components, base
import numpy as np
import time
import os

hd5file = os.environ['DATAPATH']+hd5file
io_filter = components.PyTablesSource(hd5file, dataset, f_filter=None)
base.features.Provide("SignalSource",      io_filter)
base.features.Provide("EventsOutput",      io_filter)
base.features.Provide("SpikeMarkerSource", components.SpikeDetector(contact=contact, 
                                                                    thresh=thresh,
                                                                    type=type,
                                                                    sp_win=sp_win,
                                                                    resample=10,
                                                                    align=True))
base.features.Provide("SpikeSource",       components.SpikeExtractor(sp_win=sp_win))
base.features.Provide("FeatureSource",     components.FeatureExtractor())
base.features.Provide("LabelSource",       components.ClusterAnalyzer("gmm", 4))


browser = components.SpikeBrowserWithLabels()
plot1 = components.PlotFeaturesTimeline()
plot2 = components.PlotSpikes()
legend = components.Legend()
export = components.ExportCells()

#############################################################
# Add the features here: 

base.features["FeatureSource"].add_feature("P2P")
base.features["FeatureSource"].add_feature("PCs", ncomps=2)
