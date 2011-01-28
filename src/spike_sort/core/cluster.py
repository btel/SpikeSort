#!/usr/bin/env python
#coding=utf-8

import numpy as np

def _metric_euclidean(data1, data2):
    n_pts1, n_dims1 = data1.shape
    n_pts2, n_dims2 = data2.shape
    if not n_dims1 == n_dims2:
        raise TypeError, "data1 and data2 must have the same number of columns"
    delta = np.zeros((n_pts1, n_pts2),'d')
    for d in xrange(n_dims1):
        _data1 = data1[:,d]
        _data2 = data2[:,d]
        _delta  = np.subtract.outer(_data1, _data2)**2
        delta += _delta
    return np.sqrt(delta)

def dist_euclidean(spike_waves1, spike_waves2=None):
    """Given spike_waves calculate pairwise Euclidean distance between
    them"""

    sp_data1 = spike_waves1['data']
    if spike_waves2 is None:
        sp_data2 = sp_data1
    else:
        sp_data2 = spike_waves2['data']

    d = _metric_euclidean(sp_data1.T, sp_data2.T)

    return d

def kmeans(features, n_clusters):
    """
    Automatically cluster spikes using K means algorithm
    
    :arguments:
     * features -- spike features datastructure
     * n_clusters -- number of clusters to identify
     
    :output:
     * labels
    """
    
    data = features['data']
    cl = _k_means(data, n_clusters)
    return cl

def _k_means(data, K):
    """
    Perform K means clustering
    
    :arguments:
     * data -- data vectors (n,m) where n is the number of datapoints and m is the number of variables
     * K -- number of distinct clusters to identify
     
    :output:
     * partition -- vector of cluster labels (ints) for each datapoint from `data`
    """
    
    n_dim = data.shape[1]
    centers = np.random.rand(K, n_dim)
    centers_new = np.random.rand(K, n_dim)
    partition = np.zeros(data.shape[0], dtype=np.int)
    while not (centers_new == centers).all():
        centers = centers_new.copy()
    
        distances = (centers[:,np.newaxis,:] - data)
        distances *= distances
        distances = distances.sum(axis=2)
        partition = distances.argmin(axis=0)

        for i in range(K):
            if np.sum(partition==i)>0:
                centers_new[i, :] = data[partition==i, :].mean(0)
    return partition
