from matplotlib.widgets import Lasso
import matplotlib.mlab
from matplotlib.nxutils import points_inside_poly
from matplotlib.colors import colorConverter
from matplotlib.collections import RegularPolyCollection, LineCollection

from matplotlib.pyplot import figure, show
from numpy import nonzero
from numpy.random import rand

import numpy as np
import matplotlib.pyplot as plt
import plotting

class LassoManager:
    def __init__(self, ax, data,labels, color_on='r', color_off='k'):
        self.axes = ax
        self.canvas = ax.figure.canvas
        self.data = data
        self.call_list = []

        self.Nxy = data.shape[0]
        self.color_on = colorConverter.to_rgba(color_on)
        self.color_off = colorConverter.to_rgba(color_off)

        facecolors = [self.color_on for d in range(self.Nxy)]
        fig = ax.figure
        self.collection = RegularPolyCollection(
            fig.dpi, 6, sizes=(1,),
            facecolors=facecolors,
            edgecolors=facecolors,
            offsets = self.data,
            transOffset = ax.transData)

        ax.add_collection(self.collection)

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.ind = None

    def register(self, callback_func):

        self.call_list.append(callback_func)
    def callback(self, verts):
        facecolors = self.collection.get_facecolors()
        edgecolors = self.collection.get_edgecolors()
        ind = nonzero(points_inside_poly(self.data, verts))[0]
        for i in range(self.Nxy):
            if i in ind:
                facecolors[i] = self.color_on 
                edgecolors[i] = self.color_on 
            else:
                facecolors[i] = self.color_off
                edgecolors[i] = self.color_off

        self.canvas.draw_idle()
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso
        self.ind = ind

        for func in self.call_list:
            func(ind)
    def onpress(self, event):
        if self.canvas.widgetlock.locked(): return
        if event.inaxes is None: return
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)

class SpikeWaveform:

    def __init__(self, fig, sp_dict, color_on="r", color_off="k",
                 contact=0):

        self.color_on = color_on
        self.alpha = 1.
        self.color_off = color_off
        self.figure = fig

        self.canvas = self.figure.canvas
        
        #self.lines = ax.plot(self.time, self.spikes, 'k', alpha=0.2)
        self.lines = plotting.plot_spikes(sp_dict, plot_avg=False, ec=color_off)
        #self.lines = ax.collections[0]
        
        self.Nxy = sp_dict['data'].shape[1]
        self.all_spikes = sp_dict['data']
        self.time = sp_dict['time']
        self.marked_lines = None
        
    def callback(self, ind):
        colors = np.repeat(self.color_off, self.Nxy)
        colors[ind] = self.color_on
        for lines in self.lines:
             lines.set_color(colors.tolist())
        #if self.marked_lines:
        #    self.ax.collections.remove(self.marked_lines)
        #segs = np.empty((len(ind), len(self.time), 2))
        #segs[:,:,0] = self.time[np.newaxis,:]
        #segs[:,:,1] = self.all_spikes[:,ind].T
        #self.marked_lines = LineCollection(segs,colors=self.color_on,
        #                                   alpha=self.alpha,
        #        zorder=10)
        #self.ax.add_collection(self.marked_lines)

        self.canvas.draw_idle()

def cluster_spt(spt_dict, idx):
    """return the spike times belonging to the cluster and the rest"""

    spt = spt_dict['data']
    clust = spt[idx==1]
    rest = spt[idx==0]

    return {"data": clust}, {"data":rest}


def show(features_dict, sp_dict, feat_idx,show_spikes=True):
    
    features = features_dict['data']
    names = features_dict['names']
    if type(feat_idx[0]) is int:
        ii = np.array(feat_idx)
    else:
        ii = np.array([np.nonzero(names==f)[0][0] for f in feat_idx])
        print ii
    fig_cluster = figure(figsize=(6,6))
    ax_cluster = fig_cluster.add_subplot(111, xlim=(-0.1,1.1), ylim=(-0.1,1.1),
            autoscale_on=False)
    lman = LassoManager(ax_cluster, features[:,ii], names[ii])


    if show_spikes:
        fig_spikes = figure(figsize=(6,6))
        sp_wave_plot = SpikeWaveform(fig_spikes, sp_dict)
        lman.register(sp_wave_plot.callback)

    plt.show()

    n_spikes = features.shape[0]
    clust_idx = np.zeros(n_spikes)
    clust_idx[lman.ind]=1
    return clust_idx
