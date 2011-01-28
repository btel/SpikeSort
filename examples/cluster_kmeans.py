#!/usr/bin/env python
#coding=utf-8

"""
Based on raw recordings detect spikes, calculate features and do automatic 
clustering with k-means.

After clustering the spike times are exported back to HDF5 (cell_kmeansX, where 
X is cluster index)
"""

import numpy as np
import os, sys

import spike_sort as sort
import spike_sort.io.bakerlab
import spike_sort.io.hdf5
import spike_sort.ui.manual_sort
import tables

import time

DATAPATH = "../data" 

if __name__ == "__main__":
    h5_fname = os.path.join(DATAPATH, "test_blosc.h5")
    h5f = tables.openFile(h5_fname, 'a')

    dataset = "/Gollum/s5gollum01/el3"
    sp_win = [-0.2, 0.8]
    
    start = time.time()
    sp = sort.io.hdf5.read_sp(h5f, dataset)
    spt = sort.extract.detect_spikes(sp,  contact=3,
                                     thresh=300)
    
    spt = sort.extract.align_spikes(sp, spt, sp_win, type="max", resample=10)
    sp_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    features = sort.features.combine(
            (
            sort.features.fetP2P(sp_waves),
            sort.features.fetPCs(sp_waves)),
            normalize=True
    )


    clust_idx = sort.cluster.kmeans(features,4)
    
    spike_sort.ui.plotting.plot_features(features, clust_idx)
    spike_sort.ui.plotting.figure()
    spike_sort.ui.plotting.plot_spikes(sp_waves, clust_idx,n_spikes=200)
    
    
    spike_sort.ui.plotting.show()

    #TODO: export
    #sort.io.hdf5.write_spt(clust, h5f, cell_node+"_clust",
    #                           overwrite=True)
    #sort.io.hdf5.write_spt(rest, h5f, cell_node+"_rest",
    #                           overwrite=True)

