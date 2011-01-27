#!/usr/bin/env python
#coding=utf-8
"""
Provides functions to calculate spike waveforms features.

Functions starting with `fet` implement various features calculated
from the spike waveshapes. They have usually one required argument
:ref:`spikewave` structure.  

"""


import numpy as np
import matplotlib.pyplot as plt

def combine(args, normalize=True):
    features, names = zip(*args)

    data  = np.hstack(features)
    if normalize:
        data = (data-data.min(0)[np.newaxis,:])
        data = data/data.max(0)[np.newaxis, :]

    ret_dict = {"data": data,
                "names":np.concatenate(names)}
    return ret_dict

def PCA(data,ncomps=2):
    #norm=data/np.std(data,1)[:,np.newaxis]
    #norm[np.isnan(norm)]=0
    #norm = data
    K=np.cov(data)
    evals,evecs=np.linalg.eig(K)
    order=np.argsort(evals)[::-1]
    evecs=np.real(evecs[:,order])
    evals=np.real(evals[order])
    score= np.dot(evecs[:,:ncomps].T/np.sqrt(evals[:ncomps,
                                                   np.newaxis]),data)
    return evals,evecs,score

def fetPCs(spikes_data,ncomps=2):
    """Calculate principal components (PCs)."""
   
    spikes = spikes_data["data"]
    def _getPCs(data):
        s,v,sc=PCA(data[:,:],ncomps)
        sc=(sc).astype(int)
        return sc

    if spikes.ndim==3:
        sc=[_getPCs(sp_contact.T) for sp_contact in spikes.swapaxes(0,2)]
        sc=np.vstack(sc)
        n_channels = spikes.shape[2]
    else:
        sc=_getPCs(spikes)
        n_channels = 1
    sc=sc.T
    
    names = ["Ch%d:PC%d" % (j,i) for i in range(ncomps) for j in
            range(n_channels)]
    
    return sc, names

def fetP2P(spikes_data):
    """Calculate peak-to-peak amplitudes of spike waveforms.

    :arguments:
     
     * spikes -- spikewave structure with spike waveshapes (see
       documentation for detailed specification)

    :output:

     * p2p (int) -- 2D array of peak-to-peak amplitudes in subsequent
       channels (contacts)
     
     * name -- feature labels (ex. Ch0:P2P, Ch1:P2P)

    :example:

     We will generate a spikewave structure containing only a single
     spike on a single channel

     >>> time = np.arange(0,2*np.pi,0.01) 
     >>> spikes = np.sin(time)[:,np.newaxis, np.newaxis]
     >>> spikewave = {"data": spikes, "time":time, "contacts":1, "FS":1}
     >>> p2p, labels = fetP2P(spikewave)
     >>> print p2p
     [[ 1.99999683]]

    """

    spikes = spikes_data["data"]
    p2p=spikes.max(axis=0)-spikes.min(axis=0)
    if p2p.ndim<2:
        p2p=p2p[:,np.newaxis]

    names = ["Ch%d:P2P" % i for i in range(p2p.shape[1])]

    return p2p, names

def fetSpIdx(spikes_data):
    spikes = spikes_data["data"]

    n_datapts = spikes.shape[1]

    return np.arange(n_datapts)[:, np.newaxis], ["SpIdx"]

def fetSpTime(spt_dict):

    spt = spt_dict['data']

    return spt[:, np.newaxis], ["SpTime"]
