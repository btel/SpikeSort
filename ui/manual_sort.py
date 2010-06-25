from matplotlib.widgets import Lasso
import matplotlib.mlab
from matplotlib.nxutils import points_inside_poly
from matplotlib.colors import colorConverter
from matplotlib.collections import RegularPolyCollection

from matplotlib.pyplot import figure, show
from numpy import nonzero
from numpy.random import rand

import numpy as np
import matplotlib.pyplot as plt

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
            fig.dpi, 6, sizes=(100,),
            facecolors=facecolors,
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
        ind = nonzero(points_inside_poly(self.data, verts))[0]
        for i in range(self.Nxy):
            if i in ind:
                facecolors[i] = self.color_on 
            else:
                facecolors[i] = self.color_off

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

    def __init__(self, ax, sp_dict, color_on="r", color_off="k"):

        self.color_on = colorConverter.to_rgba(color_on)
        self.color_off = colorConverter.to_rgba(color_off)

        self.canvas = ax.figure.canvas
        
        self.spikes =  sp_dict['data']
        self.time = sp_dict['time']
        self.lines = ax.plot(self.time, self.spikes, 'k', alpha=0.2)
        
        self.Nxy = self.spikes.shape[1]

    def callback(self, ind):

        for i in range(self.Nxy):
            if i in ind:
                self.lines[i].set_color(self.color_on)
            else:
                self.lines[i].set_color(self.color_off)

        self.canvas.draw_idle()

def cluster_spt(spt_dict, idx):
    """return the spike times belonging to the cluster and the rest"""

    spt = spt_dict['data']
    rest = np.asarray([spt[i] for i in range(len(spt)) if i not in idx])
    clust = spt[idx]

    return {"data": clust}, {"data":rest}


def show(features_dict, sp_dict, feat_idx):
    
    features = features_dict['data']
    names = features_dict['names']
    ii = np.array(feat_idx)
    fig = figure(figsize=(12,6))
    ax = fig.add_subplot(121, xlim=(-0.1,1.1), ylim=(-0.1,1.1),
            autoscale_on=False)
    ax2 = fig.add_subplot(122)

    sp_wave_plot = SpikeWaveform(ax2, sp_dict)
    lman = LassoManager(ax, features[:,ii], names[ii])
    lman.register(sp_wave_plot.callback)

    plt.show()

    return lman.ind
