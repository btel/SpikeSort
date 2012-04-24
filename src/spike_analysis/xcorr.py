#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

def raise_exception(*args, **kwargs):
    raise NotImplementedError("This function requires NeuroTools")

try:
    from NeuroTools.analysis import crosscorrelate
except ImportError:
    crosscorrelate = raise_exception

def show_xcorr(cells):
    n = len(cells)
    maxlag = 3
    bins = np.arange(-maxlag, maxlag, 0.01)
    ax = None
    for i in range(n):
        for j in range(i, n):
            if i != j:
                ax = plt.subplot(n, n, i + j * n + 1, sharey=ax)

                crosscorrelate(cells[i]['spt'], cells[j]['spt'], maxlag,
                                display=ax,
                                kwargs={"bins": bins})
                plt.axvline(0, color='k')

                if ax:
                    ax.set_ylabel("")
                    ax.set_yticks([])
                    ax.set_xticks([])
                    ax.set_xlabel("")
            else:
                ax_label = plt.subplot(n, n, i + j * n + 1, frameon=False)
                ax_label.set_xticks([])
                ax_label.set_yticks([])
                ax_label.text(0,0, cells[j]['dataset'],
                              transform=ax_label.transAxes,
                              rotation=45, ha='left', va='bottom')
    
    ax.set_xticks((-maxlag, maxlag))
    ax.set_yticks(ax.get_ylim())
