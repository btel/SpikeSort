#!/usr/bin/env python
#coding=utf-8

import numpy as np
from scipy import stats
import basic

import matplotlib.pyplot as plt
from io_tools import read_dataset

def plot_psth(ax, dataset, **kwargs):
    spt = dataset['spt']
    stim = dataset['stim']
    ev = dataset['ev']

    basic.plotPSTH(spt, stim,ax=ax, **kwargs)
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
    
    basic.plotraster(spt, stim,ax=ax, **kwargs)
    ymin, ymax = plt.ylim()
    if len(ev)>0:
        plt.vlines(ev, ymin, ymax)

def plot_nspikes(ax, dataset, win=[0,30], color="k"):
    spt = dataset['spt']
    stim = dataset['stim']
    
    trains = basic.SortSpikes(spt, stim, win)
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


def plot_dataset(dataset, fig=None, **kwargs):

    if not fig:
        fig = plt.gcf()
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

def show_cell(filter, cell):

    dataset = read_dataset(filter, cell)
    plot_dataset(dataset)
