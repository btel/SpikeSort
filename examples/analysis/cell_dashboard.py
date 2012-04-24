#!/usr/bin/env python
#coding=utf-8

import os
import matplotlib.pyplot as plt

from spike_analysis import dashboard
from spike_sort.io.filters import BakerlabFilter

if __name__ == "__main__":
    cell = "/Gollum/s39gollum03/el1/cell1"
    path = os.path.join(os.pardir, os.pardir, 'data', 'gollum_export.inf')
    filt = BakerlabFilter(path)
    
    dashboard.show_cell(filt, cell)

    plt.show()

    

