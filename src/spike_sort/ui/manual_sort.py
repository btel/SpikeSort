import time

from matplotlib.widgets import Lasso
from matplotlib.path import Path
from matplotlib.colors import colorConverter
from matplotlib.collections import RegularPolyCollection  # , LineCollection
from matplotlib import gridspec
from matplotlib.pyplot import figure
from numpy import nonzero

import numpy as np


class LassoManager(object):
    def __init__(self, ax, data, labels=None, color_on='r', color_off='k', markersize=1):
        self.axes = ax
        self.canvas = ax.figure.canvas
        self.data = data
        self.call_list = []

        self.Nxy = data.shape[0]
        self.color_on = colorConverter.to_rgba(color_on)
        self.color_off = colorConverter.to_rgba(color_off)

        facecolors = [self.color_on for _ in range(self.Nxy)]
        fig = ax.figure
        self.collection = RegularPolyCollection(
            fig.dpi, 6, sizes=(markersize,),
            facecolors=facecolors,
            edgecolors=facecolors,
            offsets=data,
            transOffset=ax.transData)

        ax.add_collection(self.collection, autolim=True)
        ax.autoscale_view()

        if labels is not None:
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.ind = None
        self.canvas.draw()

    def register(self, callback_func):
        self.call_list.append(callback_func)

    def callback(self, verts):
        facecolors = self.collection.get_facecolors()
        edgecolors = self.collection.get_edgecolors()
        ind = nonzero(Path(verts).contains_points(self.data))[0]
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
        if self.canvas.widgetlock.locked():
            return
        if event.inaxes is None:
            return
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata),
                           self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)


def manual_sort(features_dict, feat_idx):

    features = features_dict['data']
    names = features_dict['names']
    if type(feat_idx[0]) is int:
        ii = np.array(feat_idx)
    else:
        ii = np.array([np.nonzero(names == f)[0][0] for f in feat_idx])
    return _cluster(features[:, ii], names[:, ii])


def _cluster(data, names=None, markersize=1):
    fig_cluster = figure(figsize=(6, 6))
    gs = gridspec.GridSpec(2,2, width_ratios=[4,1], height_ratios=[4,1],
                           hspace=0, wspace=0)
    ax_cluster = fig_cluster.add_subplot(gs[0,0],
                                         xlim=(-0.1, 1.1),
                                         ylim=(-0.1, 1.1),
                                         autoscale_on=True)
    ax_marginal_x = fig_cluster.add_subplot(gs[1,0], sharex=ax_cluster,
                                            autoscale_on=False,
                                            frame_on=False)
    ax_marginal_y = fig_cluster.add_subplot(gs[0,1], sharey=ax_cluster,
                                            autoscale_on=False,
                                            frame_on=False)

    nx, binx = np.histogram(data[:, 0], 40) 
    ny, biny = np.histogram(data[:, 1], 40) 
    nx = np.append(nx, [0])
    ny = np.append(ny, [0])

    ax_marginal_x.plot(binx, nx, 'k', drawstyle='steps', clip_on=False)
    ax_marginal_x.set_ylim((nx.max(), 0))
    ax_marginal_x.set_xticks([])
    ax_marginal_x.set_yticks([])
    ax_marginal_y.plot(ny, biny, 'k', drawstyle='steps', clip_on=False)
    ax_marginal_y.set_xlim((0, ny.max()))
    ax_marginal_y.set_xticks([])
    ax_marginal_y.set_yticks([])

    lman = LassoManager(ax_cluster, data, names, markersize=markersize)

    while lman.ind is None:
        time.sleep(.01)
        fig_cluster.canvas.flush_events()

    n_spikes = data.shape[0]
    clust_idx = np.zeros(n_spikes, dtype='int16')
    clust_idx[lman.ind] = 1
    return clust_idx
