#!/usr/bin/env python
#coding=utf-8

import numpy as np

import os
import spike_sort as sort
import spike_sort.io.bakerlab
import spike_sort.io.hdf5

DATAPATH = os.environ.get("DATAPATH")

if __name__ == "__main__":
    #main
    h5_fname = os.path.join(DATAPATH, "hdf5/data_microel.h5")
    dataset = "/Joy/s3349a16/el7/cell1"
    out_dir = "/Users/bartosz/SVN/bartosz/Python/PyScripts/SpikeSorting/Data/find_missed_spikes/25-06-2010/Sim_9"
    sp_win = [-0.2, 0.8]

    spt = sort.io.bakerlab.read_spt(out_dir, "missed_rest_rest")
    sp = sort.io.hdf5.read_sp(h5_fname, dataset)
    
    spt = sort.extract.align_spikes(sp, spt, sp_win, type="min", resample=10)
    sp_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    sort.plotting.plot_spikes(sp_waves)

    features = sort.features.combine(
            (
            sort.features.fetSpIdx(sp_waves),
            sort.features.fetP2P(sp_waves),
            sort.features.fetPCs(sp_waves)))

    print features['data'].shape

    sort.plotting.plot_features(features)

    sort.plotting.show()

