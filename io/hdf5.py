#!/usr/bin/env python
#coding=utf-8

import tables
import numpy as np
import matplotlib.pyplot as plt

def read_sp(fname, dataset):

    h5f = tables.openFile(fname, 'r')

    electrode = "/".join(dataset.split('/')[:4])

    electrode_node = h5f.getNode(electrode)

    sp_raw = electrode_node.raw.read()
    FS = electrode_node.raw.attrs['sampfreq']
    
    return {"data": sp_raw, "FS": FS} 
