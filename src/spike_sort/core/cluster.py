#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

def _metric_euclidean(data1, data2):
    n_pts1, n_dims1 = data1.shape
    n_pts2, n_dims2 = data2.shape
    if n_dims1 <> n_dims2:
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
