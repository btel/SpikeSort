#!/usr/bin/env python
#coding=utf-8

"""
Simple raw data browser.

Keyboard shortcuts:

    +/- - zoom in/out
"""

import spike_sort as sort
from spike_sort.io.filters import PyTablesFilter
from spike_sort.ui import spike_browser
import os

DATAPATH = os.environ['DATAPATH']

if __name__ == "__main__":
    dataset = "/SubjectA/session01/el1"
    data_fname = os.path.join(DATAPATH, "tutorial.h5") 
    
    io_filter = PyTablesFilter(data_fname)
    sp = io_filter.read_sp(dataset) 
    spt = sort.extract.detect_spikes(sp, contact=3, thresh='auto')
    
    spike_browser.browse_data_tk(sp, spt, win=50)   
