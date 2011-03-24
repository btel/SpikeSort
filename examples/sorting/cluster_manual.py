#!/usr/bin/env python
#coding=utf-8

"""
Based on raw recordings detect spikes, calculate features and do
clustering by means of manual cluster-cutting.

After clustering the spike times are exported back to HDF5: 
    
    * selected cluster is exported to  {dataset}/cell1_clust
    * the rest is exported to  {dataset}/cell1_rest
"""

import numpy as np
import os, sys

import spike_sort as sort
from spike_sort.io.filters import PyTablesFilter, BakerlabFilter
import spike_sort.ui.manual_sort

import time


DATAPATH = "../data" 

if __name__ == "__main__":
    h5_fname = os.path.join(DATAPATH, "tutorial.h5")
    h5filter = PyTablesFilter(h5_fname, 'a')

    dataset = "/SubjectA/session01/el1"
    sp_win = [-0.2, 0.8]
    
    start = time.time()
    sp = h5filter.read_sp(dataset)
    spt = sort.extract.detect_spikes(sp,  contact=3,
                                     thresh=300)
    
    spt = sort.extract.align_spikes(sp, spt, sp_win, type="max", resample=10)
    sp_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    features = sort.features.combine(
            (
            sort.features.fetSpIdx(sp_waves),
            sort.features.fetP2P(sp_waves),
            sort.features.fetPCs(sp_waves)),
            norm=True
    )

    stop = time.time()

    clust_idx = sort.ui.manual_sort.manual_sort(features,
                                         ['Ch0:P2P','Ch3:P2P'])

    clust, rest = sort.cluster.split_cells(spt, clust_idx, [1,0])
    sort.ui.plotting.plot_spikes(sp_waves, clust_idx,n_spikes=200)
    sort.ui.plotting.show()

#    if len(clust)>0:
#        cell_node = dataset+"/cell1"
#        print "Exporting to HDF5 (file %s) as %s_{clust,rest}" % (h5_fname,
#                                                     cell_node)
#        h5filter.write_spt(clust, cell_node+"_clust",
#                               overwrite=True)
#        h5filter.write_spt(rest, cell_node+"_rest",
#                               overwrite=True)
#
#    else: 
#        print "Exiting."
    h5filter.close()
