import os

import matplotlib
matplotlib.use("TkAgg")
matplotlib.interactive(True)

from spike_beans import components, base

####################################
# Adjust these fields for your needs

# data source
hdf5filename = 'tutorial.h5'
dataset = "/SubjectA/session01/el1"

# spike detection/extraction properties
contact = 3
type = "max"
thresh = "auto"
filter_freq = (800.0, 100.0)

sp_win = [-0.6, 0.8]

path = filter(None, os.environ['DATAPATH'].split(os.sep)) + [hdf5filename]
hdf5file = os.path.join(os.sep, *path)

io_filter = components.PyTablesSource(hdf5file, dataset, f_filter=filter_freq)
base.features.Provide("SignalSource", io_filter)
base.features.Provide("EventsOutput", io_filter)
base.features.Provide("SpikeMarkerSource",
                      components.SpikeDetector(contact=contact, 
                                               thresh=thresh,
                                               type=type,
                                               sp_win=sp_win,
                                               resample=1,
                                               align=True))
base.features.Provide("SpikeSource", components.SpikeExtractor(sp_win=sp_win))
base.features.Provide("FeatureSource", components.FeatureExtractor())
base.features.Provide("LabelSource", components.ClusterAnalyzer("gmm", 4))

browser = components.SpikeBrowserWithLabels()
plot1 = components.PlotFeaturesTimeline()
plot2 = components.PlotSpikes()
legend = components.Legend()
export = components.ExportCells()

#############################################################
# Add the features here: 
base.features["FeatureSource"].add_feature("P2P")
base.features["FeatureSource"].add_feature("PCs", ncomps=2)

#############################################################
# Run the analysis (this can take a while)
browser.update()
