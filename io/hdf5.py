#!/usr/bin/env python
#coding=utf-8

import tables
import numpy as np
import matplotlib.pyplot as plt

def _get_attrs(node):
    PYTABLES_ATTRS = ["VERSION", "TITLE", "FLAVOR", "CLASS"]
    extra_attrs=dict([(name, node.getAttr(name)) 
        for name in node.attrs._v_attrnames
        if name not in PYTABLES_ATTRS])

    return extra_attrs

def read_sp(fname, dataset):

    h5f = tables.openFile(fname, 'r')

    electrode = "/".join(dataset.split('/')[:4])

    electrode_node = h5f.getNode(electrode)

    sp_raw = electrode_node.raw.read()
    FS = electrode_node.raw.attrs['sampfreq']

    h5f.close()
    
    return {"data": sp_raw, "FS": FS} 

def read_spt(fname, dataset):

    h5f = tables.openFile(fname, 'r')

    cell_node = h5f.getNode(dataset)

    spt = cell_node.read()[:].copy()
    
    ret_dict = {"data" :spt}
    extra_attrs = _get_attrs(cell_node)
    ret_dict.update(extra_attrs)
   
    cell_node.flush()
    h5f.close()
    
    return ret_dict 

def write_spt(spt_dict, fname, dataset,overwrite=False):
    
    h5f = tables.openFile(fname, 'a')
    
    spt = spt_dict['data']
    
    parts = dataset.split('/')
    group = '/'.join(parts[:-1])
    node_name = parts[-1]

    if overwrite:
        try:
            h5f.removeNode(group, node_name)
        except tables.exceptions.NodeError:
            pass

    arr_node = h5f.createArray(group, node_name, spt, 
            title="SpikeTime" , createparents=True)

    attrs = spt_dict.copy()
    del attrs['data']

    for k, v in attrs.items():
        arr_node.setAttr(k, v)

    #arr_node._v_attrs.__dict__.update(attrs)
    
    arr_node.flush()
    h5f.close()
