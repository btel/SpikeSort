#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

import os
import numpy as np
import spike_sort as sort
import spike_sort.io.hdf5
from spike_sort.ui.manual_detect import find_spikes

def calc_iso_score(dist, lam, n_spikes):
    distSS = dist[:, :n_spikes]
    distSN = dist[:, n_spikes:] 

    d0 = distSS.mean()
    expSS = np.exp(-distSS*lam*1./d0);
    expSN = np.exp(-distSN*lam*1./d0);
 
    sumSS = np.sum(expSS - np.eye(n_spikes),1)
    sumSN = np.sum(expSN, 1) 

    correctProbS = sumSS /  (sumSS + sumSN)
    isolation_score = correctProbS.mean()

    return isolation_score


if __name__ == "__main__":
    
    DATAPATH = os.environ.get("DATAPATH")
    h5_fname = os.path.join(DATAPATH, "hdf5/data_microel.h5")
    sp_win = [-0.2, 0.8]
    response_win = [8., 13.]
    
    dataset = "/Joy/s3349a16/el7/cell1_corrected"
    #dataset = "/Poppy/s32103h/el8/cell1"
    spike_type = "negative"
    stim_node = "/".join(dataset.split('/')[:3]+['stim'])
    
    sp = sort.io.hdf5.read_sp(h5_fname, dataset)
    spt = sort.io.hdf5.read_spt(h5_fname, dataset)
    stim = sort.io.hdf5.read_spt(h5_fname, stim_node)
    spt = sort.extract.align_spikes(sp, spt, sp_win, 'min')
     
    spike_waves = sort.extract.extract_spikes(sp, spt, sp_win)
    noise_waves = sort.evaluate.extract_noise_cluster(sp, spt, sp_win,
            spike_type)
    n_spikes = spike_waves['data'].shape[1]
    #i = np.random.randint(spike_waves['data'].shape[1], size=5000)
    #noise_waves = {'data':spike_waves['data'][:,i],
    #        'time':spike_waves['time']}
    all_waves = {'data': np.hstack((spike_waves['data'],
                                    noise_waves['data']))}
    dist_matrix = sort.cluster.dist_euclidean(spike_waves, all_waves)
   
    all_waves_split = spike_waves.copy()
    spike_waves_split = spike_waves.copy()
    i = np.random.rand(n_spikes).argsort()
    all_waves_split['data'] = spike_waves['data'][:, i]
    spike_waves_split['data'] = all_waves_split['data'][:, :n_spikes/2]
    
    dist_split = sort.cluster.dist_euclidean(spike_waves_split,
            all_waves_split)
    
    lambdas = np.linspace(0, 100, 20)
    iso_scores=[]
    iso_scores_split=[]
    for lam in lambdas:
         iso = calc_iso_score(dist_matrix, lam, n_spikes)
         iso_scores.append(iso)
         iso_split = calc_iso_score(dist_split, lam, n_spikes/2) 
         iso_scores_split.append(iso_split)


    plt.plot(lambdas, iso_scores, 'b')
    plt.ylabel("isolation score (spikes vs noise)")
    plt.twinx()
    plt.plot(lambdas, iso_scores_split, 'r')
    plt.ylabel("isolation score (spikes vs spikes)")
    plt.xlabel("Lambda")
    plt.title("test_iso_lambda.py")
    plt.savefig("test_iso_lambda.png")
    plt.show()
