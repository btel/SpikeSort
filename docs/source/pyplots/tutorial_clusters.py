#!/usr/bin/env python
#coding=utf-8

from spike_sort.io.filters import PyTablesFilter
from spike_sort import extract
from spike_sort import features
from spike_sort import cluster
from spike_sort.ui import plotting

dataset = '/SubjectA/session01/el1'
io_filter = PyTablesFilter('../../../data/tutorial.h5')
raw = io_filter.read_sp(dataset)
spt = extract.detect_spikes(raw,  contact=3, thresh='auto')

sp_win = [-0.2, 0.8]
spt = extract.align_spikes(raw, spt, sp_win, type="max", resample=10)
sp_waves = extract.extract_spikes(raw, spt, sp_win)
sp_feats = features.combine(
     (
      features.fetP2P(sp_waves),
      features.fetPCs(sp_waves)
     )
)
   
clust_idx = cluster.cluster("gmm",sp_feats,4)
plotting.plot_features(sp_feats, clust_idx)
plotting.show()
io_filter.close()
