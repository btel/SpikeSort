from SpikeBeans import base, components
import numpy as np
import time

dataset = "/Gollum/s39gollum01/el1"
conf_file = "../data/gollum.inf"
sp_win = [-0.2, 0.8]

base.features.Provide("SignalSource",      components.BakerlabSource(conf_file, dataset))
base.features.Provide("SpikeMarkerSource", components.SpikeDetector(contact=3, thresh='auto'))
base.features.Provide("SpikeSource",       components.SpikeExtractor(sp_win=sp_win))
base.features.Provide("FeatureSource",     components.FeatureExtractor())
base.features.Provide("LabelSource",       components.ClusterAnalyzer("gmm", 5))

base.features["FeatureSource"].add_feature("P2P")
base.features["FeatureSource"].add_feature("PCs", ncomps=1)

start = time.time()
signal = base.features['SignalSource'].signal
print time.time() - start, "s"

start = time.time()
cluster = base.features['LabelSource']
print cluster.labels[:10]
print time.time() - start, "s"
base.features['SpikeMarkerSource'].contact=2
start = time.time()
base.features['SpikeMarkerSource'].update()
print cluster.labels[:10]
print time.time() - start, "s"