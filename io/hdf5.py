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

def read_spt(fname, dataset):

    h5f = tables.openFile(fname, 'r')

    cell_node = h5f.getNode(dataset)

    spt = cell_node.read()
    
    return spt 

def wite_spt(fname, dataset, spt file=None, array_node=None, title=""):
    
    assert file
    assert array_node
    
    h5f = tables.openFile(fname, 'w')
    
    array = data['data']
    
    parts = dataset.split('/')
    group = '/'.join(parts[:-1])
    node_name = parts[-1]
    file.createArray(group, node_name, spt, title=title, createparents=True)
