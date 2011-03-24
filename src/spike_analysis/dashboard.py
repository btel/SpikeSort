#!/usr/bin/env python
#coding=utf-8

import numpy as np
from scipy import stats
import tables
import os, sys
import patterns
import hdf5tools

from NeuroTools.parameters import ParameterSet as NTParameterSet

import wx
import wx.aui
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar
import matplotlib.pyplot as plt

class Plot(wx.Panel):
    def __init__(self, parent, id = -1, dpi = None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)
        self.figure = mpl.figure.Figure(dpi=dpi, figsize=(2,2))
        self.canvas = Canvas(self, -1, self.figure)
        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        sizer.Add(self.toolbar, 0 , wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)

class PlotNotebook(wx.Panel):
    def __init__(self, parent, id = -1):
        wx.Panel.__init__(self, parent, id=id)
        self.nb = wx.aui.AuiNotebook(self)
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)

    def add(self,name="plot"):
       page = Plot(self.nb)
       self.nb.AddPage(page,name)
       return page.figure




def plot_psth(ax, dataset, **kwargs):
    spt = dataset['spt']
    stim = dataset['stim']
    ev = dataset['ev']

    patterns.plotPSTH(spt, stim,ax=ax, **kwargs)
    ymin, ymax = plt.ylim()
    if len(ev)>0:
        plt.vlines(ev, ymin, ymax)
    
    ax.text(0.95, 0.9,"total n/o spikes: %d" % (len(spt),),
            transform=ax.transAxes,
            ha='right')

def plot_isi(ax, dataset, win=[0,5], bin=0.1, color='k'):
    spt = dataset['spt']
    stim = dataset['stim']

    isi = np.diff(spt)
    intvs = np.arange(win[0], win[1], bin)
    counts, bins = np.histogram(isi, intvs)
    mode,n = stats.mode(isi)
    
    ax.set_xlim(win)
    ax.set_ylim(0, counts.max())
    ax.plot(intvs[:-1], counts, color=color, drawstyle='steps-post')
    ax.axvline(mode, color='k')
    
    ax.text(0.95, 0.9,"mode: %.2f ms" % (mode,),
             transform=ax.transAxes,
            ha='right')
    ax.set_xlabel("interval (ms)")
    ax.set_ylabel("count")

def plot_trains(ax, dataset, **kwargs):
    spt = dataset['spt']
    stim = dataset['stim']
    ev = dataset['ev']
    
    patterns.plotraster(spt, stim,ax=ax, **kwargs)
    ymin, ymax = plt.ylim()
    if len(ev)>0:
        plt.vlines(ev, ymin, ymax)

def plot_nspikes(ax, dataset, win=[0,30], color="k"):
    spt = dataset['spt']
    stim = dataset['stim']
    
    trains = patterns.SortSpikes(spt, stim, win)
    n_spks = np.array([len(t) for t in trains])
    count, bins = np.histogram(n_spks, np.arange(10))
    ax.bar(bins[:-1]-0.5, count, color=color)
    ax.set_xlim((-1,10))

    burst_frac = np.mean(n_spks>1) 
    ax.text(0.95, 0.9,"%d %% bursts" % (burst_frac*100,),
             transform=ax.transAxes,
            ha='right')
    ax.set_xlabel("no. spikes")
    ax.set_ylabel("count")


def single_cell(stim, spt, fig=None, **kwargs):

    if not fig:
        fig = plt.gcf()

    dataset = {"stim": stim['data'], "spt":spt['data'], "ev":[]}

    #dataset = hdf5tools.read_hdf5_dataset(h5f, cell)

    plt.subplots_adjust(hspace=0.3)

    ax1=fig.add_subplot(2,2,1)
    plot_psth(ax1, dataset, **kwargs)
    plt.title("PSTH")
    ax2= fig.add_subplot(2,2,2)
    plot_isi(ax2, dataset, **kwargs)
    plt.title("ISIH")
    ax3=fig.add_subplot(2,2,3)
    plot_trains(ax3,dataset)
    plt.title("raster")
    ax4= fig.add_subplot(2,2,4)
    plot_nspikes(ax4, dataset, **kwargs)
    plt.title("burst order")

def all_cells(stim, spt_cells, color_func):
    cell_ids = spt_cells.keys()

    app = wx.PySimpleApp()
    frame = wx.Frame(None,-1,'Plotter')
    plotter = PlotNotebook(frame)
    for id in cell_ids:
        fig = plotter.add('cell{0}'.format(id))
        single_cell(stim, spt_cells[id],fig=fig, color=color_func(id))
    frame.Show()
    app.MainLoop()
