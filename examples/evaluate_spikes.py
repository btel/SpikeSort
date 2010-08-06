#!/usr/bin/env python
#coding=utf-8

import os

import spike_sort as sort
import spike_sort.io.hdf5
from spike_sort.ui.manual_detect import find_spikes

if __name__ == "__main__":
    
    DATAPATH = os.environ.get("DATAPATH")
    h5_fname = os.path.join(DATAPATH, "hdf5/data_microel.h5")
    sp_win = [-0.2, 0.8]
    response_win = [8., 13.]
    
    dataset = "/Joy/s3349a16/el7/cell1"
    stim_node = "/Joy/s3349a16/stim"
    
    sp = sort.io.hdf5.read_sp(h5_fname, dataset)
    spt = sort.io.hdf5.read_spt(h5_fname, dataset)
    #stim = sort.io.hdf5.read_spt(h5_fname, stim_node)

    spt_new = sort.extract.detect_spikes(sp, -1000., 'falling')
    
    spt_aligned = sort.extract.align_spikes(sp, spt_new, sp_win, 'min')
    spt_aligned = sort.extract.remove_spikes(spt_aligned, spt, [-1,1])
    sp_waves = sort.extract.extract_spikes(sp, spt_aligned, sp_win)
    
    sort.plotting.plot_spikes(sp_waves, n_spikes=200.)


    sort.plotting.show()


