import matplotlib
matplotlib.use("TkAgg")
matplotlib.interactive(True)

from spike_beans import base, components
import numpy as np
import time

dataset = "/Gollum/s39gollum03/el1"
contact = 1
type = "min"
thresh = -200
filter_freq= (1000., 800.)
#filter_freq = None

conf_file = "../../data/gollum_export.inf"
sp_win = [-0.6, 0.8]

io_filter = components.BakerlabSource(conf_file, dataset, f_filter=filter_freq)
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
base.features.Provide("LabelSource",       components.ClusterAnalyzer("gmm", 5))


browser = components.SpikeBrowserWithLabels()
plot1 = components.PlotFeaturesTimeline()
plot2 = components.PlotSpikes()
legend = components.Legend()
export = components.ExportCells()

base.features["FeatureSource"].add_feature("P2P")
base.features["FeatureSource"].add_feature("PCs", ncomps=2)

browser.show()
plot1.show()
plot2.show()
#legend.show()
