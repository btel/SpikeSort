#!/usr/bin/env python
#coding=utf-8

import os
import itertools

import matplotlib.pyplot as plt

from spike_sort.io.filters import BakerlabFilter
from spike_analysis import io_tools, xcorr

cell_pattern = "/Gollum/s4gollum*/el*/cell*"
path = os.path.join(os.pardir, os.pardir, 'data', 'gollum.inf')
filt = BakerlabFilter(path)

if __name__ == "__main__":
    cell_nodes = io_tools.list_cells(filt, cell_pattern)
    cells = [io_tools.read_dataset(filt, node) for node in cell_nodes]
    for data, name in itertools.izip(cells, cell_nodes):
        data['dataset'] = name
    xcorr.show_xcorr(cells)
    plt.show()


