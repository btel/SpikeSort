#!/usr/bin/env python
#coding=utf-8

"""
Based on raw recordings detect spikes, calculate features and do
clustering by means of manual cluster-cutting.
"""

import os

import matplotlib
matplotlib.use("TkAgg")
matplotlib.interactive(True)

import spike_sort as sort
from spike_sort.io.filters import PyTablesFilter

DATAPATH = os.environ['DATAPATH'] 

if __name__ == "__main__":
    h5_fname = os.path.join(DATAPATH, "tutorial.h5")
    h5filter = PyTablesFilter(h5_fname, 'r')

    dataset = "/SubjectA/session01/el1"
    sp_win = [-0.2, 0.8]
    
    sp = h5filter.read_sp(dataset)
    spt = sort.extract.detect_spikes(sp, contact=3, thresh=300)
    
    spt = sort.extract.align_spikes(sp, spt, sp_win, type="max", resample=10)
    sp_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    features = sort.features.combine(
            (
            sort.features.fetSpIdx(sp_waves),
            sort.features.fetP2P(sp_waves),
            sort.features.fetPCs(sp_waves)),
            norm=True
    )
    
    clust_idx = sort.ui.manual_sort.manual_sort(features,
                                         ['Ch0:P2P', 'Ch3:P2P'])
    
    clust, rest = sort.cluster.split_cells(spt, clust_idx, [1, 0])
    
    sort.ui.plotting.figure()
    sort.ui.plotting.plot_spikes(sp_waves, clust_idx, n_spikes=200)
    
    raw_input('Press any key to exit...')
    
    h5filter.close()
    exit()
