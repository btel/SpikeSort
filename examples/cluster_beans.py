from SpikeBeans import base, components
import numpy as np
import time

dataset = "/Gollum/s39gollum02/el5"
contact = 3
conf_file = "../data/gollum_export.inf"
sp_win = [-0.3, 1]

io_filter = components.BakerlabSource(conf_file, dataset)
base.features.Provide("SignalSource",      io_filter)
base.features.Provide("EventsOutput",      io_filter)
base.features.Provide("SpikeMarkerSource", components.SpikeDetector(contact=contact, 
                                                                    thresh='auto',
                                                                    resample=10))
base.features.Provide("SpikeSource",       components.SpikeExtractor(sp_win=sp_win))
base.features.Provide("FeatureSource",     components.FeatureExtractor())
base.features.Provide("LabelSource",       components.ClusterAnalyzer("gmm", 5))


plot1 = components.PlotFeaturesTimeline()
plot2 = components.PlotSpikes()
legend = components.Legend()
export = components.ExportCells()

base.features["FeatureSource"].add_feature("P2P")
base.features["FeatureSource"].add_feature("PCs", ncomps=1)

plot1.show()
plot2.show()
legend.show()
