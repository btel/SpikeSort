#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

def PCA(data,ncomps=2):
    norm=data/std(data,1)[:,newaxis]
    K=np.cov(norm)
    evals,evecs=np.linalg.eig(K)
    order=np.argsort(evals)[::-1]
    evecs=evecs[order]
    evals=evals[order]
    score= np.matrixmultiply(evecs[:,:ncomps].T,data)
    return evals,evecs,score

def fetPCs(spikes,ncomps=2):
    """Calculate principal components (PCs).
    Input:
    - spikes - 2D array of spikes (datapoints,spikes)
    Output:
    - sc (int) - 2D array of PCs scores (spikes,score)"""
    def _getPCs(data):
        s,v,sc=PCA(data[:,:],ncomps)
        sc=(sc).astype(int)
        return sc

    if spikes.ndim==3:
        sc=[_getPCs(sp_contact.T) for sp_contact in spikes.swapaxes(0,2)]
        sc=vstack(sc)
    else:
        sc=_getPCs(spikes)
    sc=sc.T
    return sc

def fetP2P(spikes):
    """Calculate peak-to-peak amplitudes of spike waveforms.
    Input:
    - spikes - 2D array of spikes (datapoints,spikes)
    Output:
    - p2p (int) - 1D array of peak-to-peak amplitudes (spikes,score)"""
    p2p=spikes.max(axis=0)-spikes.min(axis=0)
    #p2p=p2p[:,newaxis]

    return p2p

