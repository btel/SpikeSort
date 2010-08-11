#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

import os
import numpy as np
import spike_sort as sort
import spike_sort.io.hdf5
from spike_sort.ui.manual_detect import find_spikes

from spike_sort.core.evaluate import _iso_score_dist

def get_aligned_waves(dataset, max_spikes=None):

    sp = sort.io.hdf5.read_sp(h5_fname, dataset)
    spt = sort.io.hdf5.read_spt(h5_fname, dataset)
    spt = sort.extract.align_spikes(sp, spt, sp_win, 'min')
    n_spikes =  len(spt['data'])
    if max_spikes is not None and n_spikes>max_spikes:
        i = np.random.rand(n_spikes).argsort()
        spt['data'] = spt['data'][i[:max_spikes]]
    spike_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    return spike_waves

def calc_iso_errors(dist, frac, n1, n2):

    i = np.arange(n1+n2)
    
    sh_idx1 = np.random.rand(n1).argsort()[:n1*frac]
    sh_idx2 = n1+np.random.rand(n2).argsort()[:n1*frac]

    i2 = i.copy()

    i[sh_idx1] = i2[sh_idx2]
    i[sh_idx2] = i2[sh_idx1]

    dist2 = dist[:,i][i,:]
    dist2 = dist2[:n1, :] 

    return _iso_score_dist(dist2, lam, n1)


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
    
    spike_waves1 = get_aligned_waves(dataset1, max_spikes=max_spikes)
    spike_waves2 = get_aligned_waves(dataset2, max_spikes=max_spikes)

    n_spikes1 = spike_waves1['data'].shape[1]
    n_spikes2 = spike_waves2['data'].shape[1]
    all_waves = {'data': np.hstack((spike_waves1['data'],
                                    spike_waves2['data']))}
    dist_matrix = sort.cluster.dist_euclidean(all_waves)
    iso_scores = []
    fractions = np.linspace(0, 1, 20)
    for frac in fractions:
        iso = calc_iso_errors(dist_matrix, frac, n_spikes1, n_spikes2)
        iso_scores.append(iso)


    plt.plot(fractions, iso_scores)
    plt.ylabel('isolation score')
    plt.xlabel('percent of errors')
    plt.savefig("test_iso_errors.png")
    plt.text(0.8, 0.9, "$\lambda$=%.1f" % lam)
    plt.title("test_iso_errors.py")
    plt.show()
