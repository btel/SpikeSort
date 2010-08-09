#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

import os
import numpy as np
import spike_sort as sort
import spike_sort.io.hdf5
from spike_sort.ui.manual_detect import find_spikes

def get_aligned_waves(dataset):

    sp = sort.io.hdf5.read_sp(h5_fname, dataset)
    spt = sort.io.hdf5.read_spt(h5_fname, dataset)
    spt = sort.extract.align_spikes(sp, spt, sp_win, 'min')
    spike_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    return spike_waves

def split_dataset(spike_waves, frac=0.5):
    n_spikes =  spike_waves['data'].shape[1]
    i = np.random.rand(n_spikes).argsort()

    spike_data = spike_waves['data'][:,i]
    spike_waves1 = spike_waves.copy()
    spike_waves1['data'] = spike_data[:,:int(n_spikes*frac)]
    spike_waves2 = spike_waves.copy()
    spike_waves2['data'] = spike_data[:,int(n_spikes*frac):]

    return spike_waves1, spike_waves2

if __name__ == "__main__":
    
    DATAPATH = os.environ.get("DATAPATH")
    h5_fname = os.path.join(DATAPATH, "hdf5/data_microel.h5")
    sp_win = [-0.2, 0.8]
    response_win = [8., 13.]
    lam = 10
    max_spikes=5000
    
    dataset1 = "/Joy/s3349a16/el7/cell1_corrected"
    dataset2 = "/Poppy/s32103h/el8/cell1"
    spike_type = "negative"
    
    spike_waves1 = get_aligned_waves(dataset1)
    spike_waves2 = get_aligned_waves(dataset2)
    iso_scores = np.zeros((2,2))
    splitted_a, splitted_b = split_dataset(spike_waves1)
    iso_scores[0,0] = sort.evaluate.calc_isolation_score(splitted_a,
            splitted_b, spike_type, lam=lam, max_spikes=max_spikes)
    splitted_a, splitted_b = split_dataset(spike_waves2)
    iso_scores[1,1] = sort.evaluate.calc_isolation_score(splitted_a,
            splitted_b, spike_type, lam=lam, max_spikes=max_spikes)
    iso_scores[0,1] = sort.evaluate.calc_isolation_score(spike_waves1,
            spike_waves2, spike_type, lam=lam,max_spikes=max_spikes)
    
    sort.plotting.plot_spikes(spike_waves1, n_spikes=200.)
    sort.plotting.plot_spikes(spike_waves2, n_spikes=200., ec='b')
    
    print iso_scores
