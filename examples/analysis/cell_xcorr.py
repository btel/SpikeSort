#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from spike_sort.io.filters import BakerlabFilter
from spike_analysis import io_tools, xcorr

cell_pattern = "/Gollum/s4gollum01/el*/cell*"
filter = BakerlabFilter("../../data/gollum.inf") 

if __name__ == "__main__":

    cell_nodes = io_tools.list_cells(filter, cell_pattern)
    cells = [io_tools.read_dataset(filter, node) for node in cell_nodes]
    for data, name in zip(cells, cell_nodes):
        data['dataset']=name
    xcorr.show_xcorr(cells)
    plt.show()


