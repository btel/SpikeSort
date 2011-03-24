#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

import sys

from spike_analysis import dashboard
from bakerlab import read_bakerlab_dataset
from spike_sort.io.filters import BakerlabFilter

if __name__ == "__main__":
    cell = "/Gollum/s5gollum03/el5/cell2"
    filter = BakerlabFilter("data/gollum.inf") 
    
    dashboard.show_cell(filter, cell)

    plt.show()

    

