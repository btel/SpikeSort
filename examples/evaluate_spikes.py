#!/usr/bin/env python
#coding=utf-8

import os
import numpy as np
import spike_sort as sort
import spike_sort.io.hdf5
from spike_sort.ui.manual_detect import find_spikes

if __name__ == "__main__":
    
    DATAPATH = os.environ.get("DATAPATH")
    h5_fname = os.path.join(DATAPATH, "hdf5/data_microel.h5")
    sp_win = [-0.2, 0.8]
    response_win = [8., 13.]
    
    dataset = "/Joy/s3349a16/el7/cell1_corrected"
    stim_node = "/".join(dataset.split('/')[:3]+['stim'])
    
    sp = sort.io.hdf5.read_sp(h5_fname, dataset)
    spt = sort.io.hdf5.read_spt(h5_fname, dataset)
    stim = sort.io.hdf5.read_spt(h5_fname, stim_node)
    #spt = sort.extract.detect_spikes(sp, -2000, 'falling')
    spike_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    
    snr_spike =  sort.evaluate.snr_spike(spike_waves)
    isolation_score = sort.evaluate.calc_isolation_score(sp, spt,
            sp_win, 'negative')
    
    noise_waves = sort.evaluate.extract_noise_cluster(sp, spt, sp_win,
            'negative') 
    sort.plotting.plot_spikes(noise_waves, n_spikes=200.)
    sort.plotting.plot_spikes(spike_waves, n_spikes=200., ec='b')

     
    print "SNR_{spk} (spike) =", sort.evaluate.snr_spike(spike_waves)
    print "SNR_{spk} (noise) =", sort.evaluate.snr_spike(noise_waves)
    print "Isolation Score = ", isolation_score


    sort.plotting.show()


