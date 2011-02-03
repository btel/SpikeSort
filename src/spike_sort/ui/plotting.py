#!/usr/bin/env python
#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.pyplot import show, figure, close
from matplotlib.collections import LineCollection

cmap = plt.cm.jet

def label_color(labels):
    """Map labels to number range [0,1]"""
    
    num_labels = np.linspace(0,1., len(labels))
    mapper = dict(zip(labels, num_labels))
    
    @np.vectorize
    def map_func(lab):
        return mapper[lab]
    
    def color_func(lab):
        return cmap(map_func(lab))
    return color_func


def plot_spikes(spike_data, clust_idx=None,ec=None, alpha=0.2,
        n_spikes='all', contacts='all', plot_avg=True):
    """Plot Spike waveshapes

    :arguments:
    
     * spike_data : dict
     * clust_idx : sequence
       sequence of the length equal to the number of spikes; labels of
       clusters to which spikes belong
     * n_spikes : int or "all"
       number of spikes to plot; 'all' if all
     * plot_avg: bool
       plot waveform averages?

    :output:
    
     * lines_segments
       matplotlib line collection of spike waveshapes
    """

    spikes =  spike_data['data']
    time = spike_data['time']

    if contacts == 'all':
        contacts = np.arange(spikes.shape[2])
    
    n_pts = len(time)
    
    if clust_idx is None:
        clust_idx = np.ones(spikes.shape[1])
        
    norm = label_color(np.unique(clust_idx))
    unique_labs = np.unique(clust_idx)
    
    if  not n_spikes=='all':
        spikes_lab = [spikes[:,clust_idx==c,:][:,:n_spikes,:] 
                      for c in unique_labs]
        spikes = np.concatenate(spikes_lab,1)
        new_idx = [np.repeat(c, np.min([n_spikes, np.sum(c==clust_idx)])) 
                   for c in unique_labs]
        clust_idx = np.concatenate(new_idx)
    n_spikes = spikes.shape[1]

    if ec:
        colors = ec
    else:
        colors = norm(clust_idx)
    
    line_segments = []
    for i, contact_id in enumerate(contacts): 
        ax = plt.subplot(2,2, i+1)
        ax.set_xlim(time.min(), time.max())
        ax.set_ylim(spikes.min(), spikes.max())
        
        segs = np.zeros((n_spikes, n_pts, 2))
        segs[:,:,0] = time[np.newaxis,:]
        segs[:,:,1] = spikes[:,:, contact_id].T
        collection = LineCollection(segs,colors=colors,
                                            alpha=alpha)
        line_segments.append(collection)
        ax.add_collection(collection)
        
        if plot_avg:
            
            if ec is None:
                colors_mean = norm(unique_labs)
            spikes_mean = np.array([spikes[:, clust_idx==l, contact_id].mean(1) for l in
                unique_labs])
            segs_mean = np.zeros((len(unique_labs), n_pts, 2))
            segs_mean[:,:,0] = time[np.newaxis,:]
            segs_mean[:,:,1] = spikes_mean
            ls_mean = LineCollection(segs_mean, colors=colors_mean,
                    linewidths=2, zorder=10, alpha=1)
            ax.add_collection(ls_mean)
    return line_segments

def plot_features(features_dict, clust_idx=None, size=1):
    """Plot features and their histograms
    
    :arguments:
     * features_dict --features data structure (dicitionary with
       `data` and `names` keys)
     * clust_idx (default: None) -- array of size (n_spikes,)
       containing indices of clusters to which each spike was
       classfied
     * size (default: 1)
    
    """

    features = features_dict['data']
    names = features_dict['names']

    n_spikes, n_feats = features.shape
    
    if clust_idx is None:
        clust_idx = np.ones(n_spikes)
    
    norm = label_color(np.unique(clust_idx))
    fig = plt.gcf()
    for i in range(n_feats):
        for j in range(n_feats):
            ax = fig.add_subplot(n_feats, n_feats, j*n_feats + i + 1)
            if i<>j:
                for c in np.unique(clust_idx):
                    plt.plot(features[clust_idx==c,i],
                             features[clust_idx==c,j],".", 
                            color=norm(c), markersize=size)
            else:
                ax.set_frame_on(False)
                for c in np.unique(clust_idx):
                    plt.hist(features[clust_idx==c,i],20, 
                            [0,1], ec="none", fc=norm(c),
                            alpha=0.5, normed=True)
            plt.xticks([])
            plt.yticks([])
    
    for i in range(n_feats):
        ax = plt.subplot(n_feats, n_feats, i+1)
        ax.set_xlabel(names[i])
        ax.xaxis.set_label_position("top")
        ax = plt.subplot(n_feats, n_feats, i*n_feats + 1)
        ax.set_ylabel(names[i])




