#!/usr/bin/env python
#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.pyplot import show, figure, close
from matplotlib.collections import LineCollection
import spike_sort
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

def plot_spikes(spikes, clust_idx=None, show_cells='all', **kwargs):
    """Plot Spike waveshapes

    Parameters
    ----------
    spike_data : dict
    clust_idx : sequence
        sequence of the length equal to the number of spikes; labels of
        clusters to which spikes belong
    show_cells : list or 'all'
        list of identifiers of clusters (cells) to plot
    plot_avg: bool
        if True plot waveform averages

    Returns
    -------
    lines_segments : object
        matplotlib line collection of spike waveshapes
    """
    
    if clust_idx is None:
        spikegraph(spikes,**kwargs)
    else:
        spikes_cell = spike_sort.extract.split_cells(spikes, clust_idx)
        
        if show_cells == 'all':
            labs = spikes_cell.keys()
        else:
            labs = show_cells
       
        color_func = label_color(spikes_cell.keys())
        for l in labs:
            spikegraph(spikes_cell[l], color_func(l), **kwargs)
    
def spikegraph(spike_data, color='k', alpha=0.2, n_spikes='all', contacts='all', 
                plot_avg=True, fig=None):


    spikes =  spike_data['data']
    time = spike_data['time']

    if contacts == 'all':
        contacts = np.arange(spikes.shape[2])
    
    n_pts = len(time)
    
    if  not n_spikes=='all':
        spikes = spikes[:,:n_spikes,:]
    n_spikes = spikes.shape[1]
    if fig is None:
        fig = plt.gcf()
    line_segments = []
    for i, contact_id in enumerate(contacts): 
        ax = fig.add_subplot(2,2, i+1)
        #ax.set_xlim(time.min(), time.max())
        #ax.set_ylim(spikes.min(), spikes.max())
        
        
        segs = np.zeros((n_spikes, n_pts, 2))
        segs[:,:,0] = time[np.newaxis,:]
        segs[:,:,1] = spikes[:,:, contact_id].T
        collection = LineCollection(segs,colors=color,
                                            alpha=alpha)
        line_segments.append(collection)
        ax.add_collection(collection, autolim=True)
        
        if plot_avg:
            spikes_mean = spikes[:, :, i].mean(1) 
            ax.plot(time, spikes_mean, color='w',lw=3)
            ax.plot(time, spikes_mean, color=color,lw=2)
        ax.autoscale_view(tight=True)
            
    return line_segments

def plot_features(features, clust_idx=None, show_cells='all', **kwargs):
    """Plot features and their histograms
    
    Parameters
    ----------
    features_dict : dict
        features data structure
    clust_idx : array or None
        array of size (n_spikes,) containing indices of clusters to which
        each spike was classfied
    show_cells : list or 'all'
        list of identifiers of clusters (cells) to plot
    
    """

    if clust_idx is None:
        featuresgraph(features,**kwargs)
    else:
        features_cell = spike_sort.features.split_cells(features, clust_idx)
        
        if show_cells == 'all':
            labs = features_cell.keys()
        else:
            labs = show_cells
       
        color_func = label_color(features_cell.keys())
        for l in labs:
            featuresgraph(features_cell[l], color_func(l),**kwargs)
    
def featuresgraph(features_dict, color='k', size=1, datarange=None, fig=None):

    features = features_dict['data']
    names = features_dict['names']

    _, n_feats = features.shape
    if fig is None:
        fig = plt.gcf()
    axes = [[fig.add_subplot(n_feats, n_feats, i*n_feats + j + 1) 
             for i in range(n_feats)] for j in range(n_feats)]
    for i in range(n_feats):
        for j in range(n_feats):
            ax = axes[i][j]
            if i!=j:
                ax.plot(features[:,i],
                         features[:,j],".", 
                        color=color, markersize=size)
                if datarange:
                    ax.set_xlim(datarange)
                    
            else:
                ax.set_frame_on(False)
           
                ax.hist(features[:,i],20, 
                        datarange, ec="none", fc=color,
                        alpha=0.5, normed=True)
                if datarange:
                    ax.set_xlim(datarange)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel(names[i])
            ax.set_ylabel(names[j])
            ax.xaxis.set_label_position("top")
            ax.xaxis.label.set_visible(False)
            ax.yaxis.label.set_visible(False)

    for i in range(n_feats):
        ax = axes[i][0] 
        ax.xaxis.label.set_visible(True)
        ax = axes[0][i]
        ax.yaxis.label.set_visible(True)

def legend(labels, colors=None, ax=None):
    
    if ax is None:
        ax=plt.gca()
    if colors is None:
        color_func = label_color(labels)
        colors = [color_func(i) for i in labels]
        
    ax.set_frame_on(False)
    n_classes = len(labels)
    x, y =np.ones(n_classes), np.arange(n_classes)
    ax.scatter(x, y, c=colors, marker='s',edgecolors="none",s=100)
    ax.set_xticks([])
    ax.set_yticks([])
    for i, l in enumerate(labels):
        ax.text(x[i]+0.005, y[i], "Cell {0}".format(l), va='center', ha='left', 
                 transform=ax.transData)
    

