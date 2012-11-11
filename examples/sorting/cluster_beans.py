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
detection_type = "max"
thresh = "auto"
filter_freq = (800.0, 100.0)

sp_win = [-0.6, 0.8]

path = filter(None, os.environ['DATAPATH'].split(os.sep)) + [hdf5filename]
hdf5file = os.path.join(os.sep, *path)

io = components.PyTablesSource(hdf5file, dataset)
io_filter = components.FilterStack()

base.register("RawSource", io)
base.register("EventsOutput", io)
base.register("SignalSource", io_filter)
base.register("SpikeMarkerSource",
                      components.SpikeDetector(contact=contact, 
                                               thresh=thresh,
                                               type=detection_type,
                                               sp_win=sp_win,
                                               resample=1,
                                               align=True))
base.register("SpikeSource", components.SpikeExtractor(sp_win=sp_win))
base.register("FeatureSource", components.FeatureExtractor())
base.register("LabelSource", components.ClusterAnalyzer("gmm", 4))

browser = components.SpikeBrowserWithLabels()
feature_plot = components.PlotFeaturesTimeline()
wave_plot = components.PlotSpikes()
legend = components.Legend()
export = components.ExportCells()

#############################################################
# Add filters here:
base.features["SignalSource"].add_filter("LinearIIR", *filter_freq)

# Add the features here: 
base.features["FeatureSource"].add_feature("P2P")
base.features["FeatureSource"].add_feature("PCs", ncomps=2)

#############################################################
# Run the analysis (this can take a while)
browser.update()
