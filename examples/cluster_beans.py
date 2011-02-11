from SpikeBeans import base, components
import numpy as np
import time

dataset = "/Gollum/s39gollum01/el1"
conf_file = "../data/gollum.inf"
sp_win = [-0.2, 0.8]

class PrintLabel(base.Component):
    cluster = base.RequiredFeature("LabelSource")
    def update(self):
        start = time.time()
        labels =  self.cluster.labels
        duration = time.time() - start
        
        print "Cell labels:", labels, "No. labels:", len(labels), "Calc. time:", duration
        

base.features.Provide("SignalSource",      components.BakerlabSource(conf_file, dataset))
base.features.Provide("SpikeMarkerSource", components.SpikeDetector(contact=3, thresh='auto'))
base.features.Provide("SpikeSource",       components.SpikeExtractor(sp_win=sp_win))
base.features.Provide("FeatureSource",     components.FeatureExtractor())
base.features.Provide("LabelSource",       components.ClusterAnalyzer("gmm", 5))
base.features.Provide("Print",             PrintLabel())

plot = components.PlotFeatures()
base.features["FeatureSource"].add_feature("P2P")
base.features["FeatureSource"].add_feature("PCs", ncomps=1)

start = time.time()
signal = base.features['SignalSource'].signal
print time.time() - start, "s"
plot.draw()

cluster = base.features['Print'].update()

base.features['SpikeMarkerSource'].contact=2
base.features['SpikeMarkerSource'].update()
plot.draw(True)
plot.show()