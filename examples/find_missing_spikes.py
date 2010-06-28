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
    trigger = "/Joy/s3349a16/stim"
    
    sp = sort.io.hdf5.read_sp(h5_fname, dataset)
    spt = sort.io.hdf5.read_spt(h5_fname, dataset)
    trigger = sort.io.hdf5.read_spt(h5_fname, trigger)

    missed_spikes = find_spikes(sp, spt, trigger, response_win, sp_win)

    sp_waves = sort.extract.extract_spikes(sp, missed_spikes, sp_win)
    sort.plotting.plot_spikes(sp_waves)

    sort.plotting.show()


