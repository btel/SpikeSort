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
import spike_sort.io.bakerlab
import spike_sort.io.hdf5
import spike_sort.ui.manual_sort
import tables

DATAPATH = "../data" 

if __name__ == "__main__":
    h5_fname = os.path.join(DATAPATH, "sample.h5")
    h5f = tables.openFile(h5_fname, 'a')

    dataset = "/Gollum/s5gollum01/el3"
    sp_win = [-0.2, 0.8]

    sp = sort.io.hdf5.read_sp(h5f, dataset)
    spt = sort.extract.detect_spikes(sp,  contact=3,
                                     thresh=300)
    
    spt = sort.extract.align_spikes(sp, spt, sp_win, type="max", resample=10)
    sp_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    features = sort.features.combine(
            (
            sort.features.fetSpIdx(sp_waves),
            sort.features.fetP2P(sp_waves),
            sort.features.fetPCs(sp_waves)),
            normalize=True)


    clust_idx = sort.ui.manual_sort.show(features, sp_waves,
                                         ['Ch0:P2P','Ch3:P2P'],
                                         show_spikes=True)

    clust, rest = sort.ui.manual_sort.cluster_spt(spt, clust_idx)

    if len(clust)>0:
        cell_node = dataset+"/cell1"
        print "Exporting to HDF5 (file %s) as %s_{clust,rest}" % (h5_fname,
                                                     cell_node)
        sort.io.hdf5.write_spt(clust, h5f, cell_node+"_clust",
                               overwrite=True)
        sort.io.hdf5.write_spt(rest, h5f, cell_node+"_rest",
                               overwrite=True)

    else: 
        print "Exiting."
