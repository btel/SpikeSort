#!/usr/bin/env python
#coding=utf-8

"""
Based on raw recordings detect spikes, calculate features and do automatic 
clustering with k-means.

"""

import os, sys
sys.path.append("../modules")
import dashboard 
import numpy as np

import spike_sort as sort
import spike_sort.ui.manual_sort
import tables

import time
from spike_sort.io.filters import PyTablesFilter, BakerlabFilter
from spike_sort.io import export

DATAPATH = "../data" 


if __name__ == "__main__":

    dataset = "/Gollum/s39gollum01/el1"
    cell_templ = "/Gollum/s39gollum01/el1/cell{cell_id}"
    stim_node = "/".join(dataset.split("/")[:3]+["stim"])

    sp_win = [-0.2, 0.8]

    start = time.time()
    io_filter = BakerlabFilter("../data/gollum.inf")
    export_filter = PyTablesFilter("../data/exported.h5")

    sp = io_filter.read_sp(dataset,memmap="numpy")
    spt = sort.extract.detect_spikes(sp,  contact=3, thresh='auto')
    
    spt = sort.extract.align_spikes(sp, spt, sp_win, type="max",
                                    contact=3, resample=10)
    sp_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    features = sort.features.combine(
            (
            sort.features.fetP2P(sp_waves,contacts=[0,1,2,3]),
            sort.features.fetPCs(sp_waves)),
            normalize=True
    )


    clust_idx = sort.cluster.cluster("gmm",features,5)
    
    features = sort.features.combine(
            (sort.features.fetSpIdx(sp_waves), features))
    spike_sort.ui.plotting.plot_features(features, clust_idx)
    spike_sort.ui.plotting.figure()
    spike_sort.ui.plotting.plot_spikes(sp_waves, clust_idx, n_spikes=200)
    
    spt_cells = sort.cluster.split_cells(spt, clust_idx)
    features_cells = sort.features.split_cells(features, clust_idx)
    spikes_cells = sort.extract.split_cells(sp_waves, clust_idx)
    stim = io_filter.read_spt(dataset)
    io_filter.close()
   
    from matplotlib.pyplot import figure,show

    color_map = spike_sort.ui.plotting.label_color(np.unique(spt_cells.keys()))
    for i in spt_cells.keys():
        #plotPSTH(spt_cells[i]['data'], stim['data'], 
        #         color=color_map(i),
        #         label="cell {0}".format(i)) 
        figure()
        dashboard.single_cell(stim, spt_cells[i], color=color_map(i))
    #figure()
    #spike_sort.ui.plotting.legend(spt_cells.keys(),
    #                             color_map(spt_cells.keys()))
    show()
    #dashboard.all_cells(stim, spt_cells, color_map)
    #export.export_cells(export_filter, cell_templ, spt)
    

