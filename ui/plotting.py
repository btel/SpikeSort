#!/usr/bin/env python
#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.pyplot import show

def plot_spikes(spike_data):

    spikes =  spike_data['data']
    time = spike_data['time']

    plt.plot(time, spikes, 'k', alpha=0.2)

def plot_features(features_dict):

    features = features_dict['data']
    names = features_dict['names']

    n_feats = features.shape[1]

    fig = plt.figure()
    for i in range(n_feats):
        for j in range(n_feats):
            ax = fig.add_subplot(n_feats, n_feats, j*n_feats + i + 1)
            if i<>j:
                plt.plot(features[:,i], features[:, j], '.')
                plt.xlabel(names[i])
                plt.ylabel(names[j])
                plt.xticks([])
                plt.yticks([])
            else:
                ax.set_frame_on(False)
                plt.hist(features[:,i])


