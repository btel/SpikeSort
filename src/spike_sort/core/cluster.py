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
    from scikits.learn import cluster as skcluster
    from scikits.learn import mixture

    def k_means_plus(*args, **kwargs):
        """k means with smart initialization.

        Notes
        -----
        This function requires scikits-learn

        See Also
        --------
        kmeans
        """
        return skcluster.k_means(*args, **kwargs)[1]

    def gmm(data, k):
        """Cluster based on gaussian mixture models

        Parameters
        ----------
        data : dict
            features structure
        k :  int
            number of clusters

        Returns
        -------
        cl : int array
            cluster indicies

        Notes
        -----
        This function requires scikits-learn
        """
        try:
            clf = mixture.GMM(n_states=k, cvtype='full')
        except TypeError:
            clf = mixture.GMM(n_components=k, cvtype='full')
        clf.fit(data)
        cl = clf.predict(data)
        return cl


except ImportError:
    k_means_plus = default_scikits("k_means_plus")
    gmm = default_scikits("gmm")


from spike_sort.ui import manual_sort


def manual(data, *args, **kwargs):
    """Sort spikes manually by cluster cutting

    Opens a new window in which you can draw cluster of arbitrary
    shape.

    Notes
    -----
    Only two first features are plotted
    """
    return manual_sort._cluster(data[:, :2])


def none(data):
    """Do nothing"""
    return np.zeros(data.shape[0], dtype='int16')


def _metric_euclidean(data1, data2):
    n_pts1, n_dims1 = data1.shape
    n_pts2, n_dims2 = data2.shape
    if not n_dims1 == n_dims2:
        raise TypeError("data1 and data2 must have the same number of columns")
    delta = np.zeros((n_pts1, n_pts2), 'd')
    for d in xrange(n_dims1):
        _data1 = data1[:, d]
        _data2 = data2[:, d]
        _delta = np.subtract.outer(_data1, _data2) ** 2
        delta += _delta
    return np.sqrt(delta)


def dist_euclidean(spike_waves1, spike_waves2=None):
    """Given spike_waves calculate pairwise Euclidean distance between
    them"""

    sp_data1 = np.concatenate(spike_waves1['data'], 1)

    if spike_waves2 is None:
        sp_data2 = sp_data1
    else:
        sp_data2 = np.concatenate(spike_waves2['data'], 1)
    d = _metric_euclidean(sp_data1, sp_data2)

    return d


def cluster(method, features, *args, **kwargs):
    """Automatically cluster spikes using K means algorithm

    Parameters
    ----------
    features : dict
        spike features datastructure
    n_clusters : int
        number of clusters to identify
    args, kwargs :
        optional arguments that are passed to the clustering algorithm

    Returns
    -------
    labels : array
        array of cluster (unit) label - one for each cell

    Examples
    --------
    Create a sample feature dataset and use k-means clustering to find
    groups of spikes (units)

    >>> import spike_sort
    >>> import numpy as np
    >>> np.random.seed(1234) #k_means uses random initialization
    >>> features = {'data':np.array([[0.,0.],
    ...                              [0, 1.],
    ...                              [0, 0.9],
    ...                              [0.1,0]])}
    >>> labels = spike_sort.cluster.cluster('k_means', features, 2)
    >>> print labels
    [0 1 1 0]
    """
    try:
        cluster_func = eval(method)
    except NameError:
        raise NotImplementedError(
                    "clustering method %s is not implemented" % method
                    )

    data = features['data']
    mask = features.get('is_valid')
    if mask is not None:
        valid_data = data[mask, :]
        cl = cluster_func(valid_data, *args, **kwargs)
        labels = np.zeros(data.shape[0], dtype='int') - 1
        labels[mask] = cl
    else:
        labels = cluster_func(data, *args, **kwargs)
    return labels


def k_means(features, K):
    """Perform K means clustering

    Parameters
    ----------
    data : dict
        data vectors (n,m) where n is the number of datapoints and m is
        the number of variables
    K : int
        number of distinct clusters to identify

    Returns
    -------
    partition : array
        vector of cluster labels (ints) for each datapoint from `data`
    """
    n_dim = features.shape[1]
    centers = np.random.rand(K, n_dim)
    centers_new = np.random.rand(K, n_dim)
    partition = np.zeros(features.shape[0], dtype=np.int)
    while not (centers_new == centers).all():
        centers = centers_new.copy()

        distances = (centers[:, np.newaxis, :] - features)
        distances *= distances
        distances = distances.sum(axis=2)
        partition = distances.argmin(axis=0)

        for i in range(K):
            if np.sum(partition == i) > 0:
                centers_new[i, :] = features[partition == i, :].mean(0)
    return partition


def split_cells(spt_dict, idx, which='all'):
    """return the spike times belonging to the cluster and the rest"""

    if which == 'all':
        classes = np.unique(idx)
    else:
        classes = which
    spt = spt_dict['data']
    spt_dicts = dict([(cl, {'data': spt[idx == cl]}) for cl in classes])
    return spt_dicts
