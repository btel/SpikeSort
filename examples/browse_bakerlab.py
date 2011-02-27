#!/usr/bin/env python
#coding=utf-8

"""
Based on raw recordings detect spikes, calculate features and do automatic 
clustering with k-means.

TODO:
After clustering the spike times are exported back to HDF5 (cell_kmeansX, where 
X is cluster index)
"""

import numpy as np
import os, sys

import spike_sort as sort
import spike_sort.ui.manual_sort
import tables

import time
from spike_sort.io.filters import PyTablesFilter, BakerlabFilter
from spike_sort.ui import spike_browser


DATAPATH = "../data" 

if __name__ == "__main__":

    dataset = "/Gollum/s39gollum02/el2"
    sp_win = [-0.2, 0.8]
    
    start = time.time()
    io_filter = BakerlabFilter("../data/gollum.inf")
    io_spt_filter = BakerlabFilter("../data/gollum_export.inf")
    filter = sort.extract.Filter('ellip', 300., 100.)
    #filter = None
    sp = sort.extract.filter_proxy(io_filter.read_sp(dataset,memmap="none"), filter)
    spt = io_spt_filter.read_spt(dataset+'/cell1')
    #spt = sort.extract.detect_spikes(sp,  contact=1,
    #                                 thresh='auto')
    
    spike_browser.browse_data(sp, spt, win=50)
    sort.plotting.show()
    
    