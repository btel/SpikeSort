#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from NeuroTools.analysis import crosscorrelate
from spike_sort.io.filters import BakerlabFilter
from spike_analysis import io_tools

cell_pattern = "/Gollum/s4gollum01/el*/cell*"
filter = BakerlabFilter("data/gollum.inf") 

if __name__ == "__main__":

    cell_nodes = io_tools.list_cells(filter, cell_pattern)

    cells = [io_tools.read_dataset(filter, node) for node in cell_nodes]
   
    n = len(cells)
    maxlag = 3
    binsz = 0.1
    bins = np.arange(-maxlag, maxlag, 0.01)
    ax=None
    for i in range(n):
        for j in range(i,n):
            if i<>j:
                ax = plt.subplot(n,n,i+j*n+1, sharey=ax)

                crosscorrelate(cells[i]['spt'], cells[j]['spt'], maxlag,
                                display=ax,
                                kwargs={"bins":bins})
                plt.axvline(0, color='k')

                if ax:
                    ax.set_ylabel("")
                    ax.set_yticks([])
                    ax.set_xticks([])
                    ax.set_xlabel("")
            else:
                ax_label = plt.subplot(n,n,i+j*n+1, frameon=False)
                ax_label.set_xticks([])
                ax_label.set_yticks([])
                ax_label.text(0,0, cell_nodes[j],
                              transform=ax_label.transAxes,
                              rotation=45, ha='left',
                              va='bottom')

        
    ax.set_xticks((-maxlag, maxlag))
    ymin, ymax = ax.get_ylim()
    ax.set_yticks((ymin, ymax))
    plt.show()


