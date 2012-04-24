#!/usr/bin/env python
#coding=utf-8

"""
Based on raw recordings detect spikes, calculate features and do automatic 
clustering with gaussian mixture models.
"""

import os

import spike_sort as sort
from spike_sort.io.filters import PyTablesFilter
import spike_sort.ui.manual_sort

DATAPATH = os.environ['DATAPATH'] 

if __name__ == "__main__":
    h5_fname = os.path.join(DATAPATH, "tutorial.h5")
    h5filter = PyTablesFilter(h5_fname, 'a')

    dataset = "/SubjectA/session01/el1"
    sp_win = [-0.2, 0.8]
    
    sp = h5filter.read_sp(dataset)
    spt = sort.extract.detect_spikes(sp, contact=3, thresh='auto')
    
    spt = sort.extract.align_spikes(sp, spt, sp_win, type="max", resample=10)
    sp_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    features = sort.features.combine(
            (sort.features.fetP2P(sp_waves),
             sort.features.fetPCs(sp_waves)),
            norm=True
    )


    clust_idx = sort.cluster.cluster("gmm",features,4)
    
    spike_sort.ui.plotting.plot_features(features, clust_idx)
    spike_sort.ui.plotting.figure()
    spike_sort.ui.plotting.plot_spikes(sp_waves, clust_idx,n_spikes=200)
    
    
    spike_sort.ui.plotting.show()
    h5filter.close()
