#!/usr/bin/env python
#coding=utf-8

import numpy as np

def default_scikits(method):
    def foo(*args):
        raise NotImplementedError(
                    "scikits.learn must be installed to use %s" % method
                    )
    return foo

try:
    import scikits.learn.cluster
    from scikits.learn import mixture
    
    def k_means_plus(*args, **kwargs):
        """k means with smart initialization.
        
        see also kmeans"""
        return scikits.learn.cluster.k_means(*args, **kwargs)[1]
    
    def gmm(data, k):
        """cluster based on gaussian mixture models (from scikits.learn)
        
        :arguments:
         * data -- features structures
         * k -- number of clusters
         
         :output:
          * k -- number of clusters to fit
        """
        clf = mixture.GMM(n_states=k, cvtype='full')
        clf.fit(data)
        cl = clf.predict(data)
        return cl

    
except ImportError:
    k_means_plus = default_scikits("k_means_plus")



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

def cluster(method, features,  *args, **kwargs):
    """
    Automatically cluster spikes using K means algorithm
    
    :arguments:
     * features -- spike features datastructure
     * n_clusters -- number of clusters to identify
     
    :output:
     * labels
    """
    try:
        
        cluster_func = eval(method)
    except NameError:
        raise NotImplementedError(
                    "clustering method %s is not implemented" % method
                    )
    
    data = features['data']
    cl = cluster_func(data, *args, **kwargs)
    return cl

def k_means(data, K):
    """
    Perform K means clustering
    
    :arguments:
     * data -- data vectors (n,m) where n is the number of datapoints and m is 
       the number of variables
     * K -- (required) number of distinct clusters to identify
     
    :output:
     * partition -- vector of cluster labels (ints) for each datapoint from 
       `data`
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

def cluster2spt(spt_dict, idx, which='all'):
    """return the spike times belonging to the cluster and the rest"""

    if which == 'all':
        classes = np.unique(idx)
    else:
        classes = which
    
    spt = spt_dict['data']
    spt_dicts = [{'data': spt[idx==cl], 'cell_id': cl} for cl in classes]

    return spt_dicts