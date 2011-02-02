#!/usr/bin/env python
#coding=utf-8

"""
HDF5 is a hierarchical datafile -- data is organised in a tree. The
standard layout is::
   
    /{SubjectName}/
    /{SubjectName}/{SessionName}/{ElectrodeID}/
    /{SubjectName}/{SessionName}/{ElectrodeID}/stim: stimulus time
    /{SubjectName}/{SessionName}/{ElectrodeID}/raw: spike waveforms
    /{SubjectName}/{SessionName}/{ElectrodeID}/{CellID}: spike waveforms
    /{SubjectName}/{SessionName}/{ElectrodeID}/{CellID}/spt: spike
    times

where curly brackets `{}` denote a group.

This layout may be adjusted by changing paths

TODO: implement custom layouts
"""

_open_files = []

def _get_attrs(node):
    PYTABLES_ATTRS = ["VERSION", "TITLE", "FLAVOR", "CLASS"]
    extra_attrs=dict([(name, node.getAttr(name)) 
        for name in node.attrs._v_attrnames
        if name not in PYTABLES_ATTRS])

    return extra_attrs

def _open_file(fname, mode='r'):
    if type(fname) is tables.File:
        h5f = fname
    else:
        try:
            h5f = tables.openFile(fname, mode)
        except ValueError:
            fname.close()
            _open_files.remove(fname)
            h5f = tables.openFile(fname, mode)
    if not h5f in _open_files:
        _open_files.append(h5f)
    return h5f

def close_all():
    global _open_files
    for f in _open_files:
        f.close()
    _open_files=[]

def read_sp(fname, dataset):
    """Read continous waveforms (EEG, LFG, spike waveform.
    
    :arguments:
        * fname -- filename or open hdf5 file object
        * dataset -- (string) path pointing to cell node
    """
    
    h5f = _open_file(fname, 'r+')
    electrode = "/".join(dataset.split('/')[:4])

    electrode_node = h5f.getNode(electrode)

    sp_raw = electrode_node.raw
    FS = electrode_node.raw.attrs['sampfreq']

    try:
        n_contacts = sp_raw.shape[1]
    except IndexError:
        n_contacts = 1
    
    return {"data": sp_raw, "FS": FS, "n_contacts":n_contacts} 

def read_spt(fname, dataset):
    """Read event times (such as spike or stimulus times).
   
    :arguments:
        * fname -- filename or open hdf5 file object
        * dataset -- (string) with path pointing to cell node
    """

    h5f = _open_file(fname, 'r+')

    cell_node = h5f.getNode(dataset)

    spt = cell_node.read()[:].copy()
    
    ret_dict = {"data" :spt}
    extra_attrs = _get_attrs(cell_node)
    ret_dict.update(extra_attrs)
   
    cell_node.flush()
    #h5f.close()
    
    return ret_dict 

def write_spt(spt_dict, fname, dataset,overwrite=False):
    
    h5f = _open_file(fname, 'a')
    
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
            title="" , createparents=True)

    attrs = spt_dict.copy()
    del attrs['data']

    for k, v in attrs.items():
        arr_node.setAttr(k, v)

    #arr_node._v_attrs.__dict__.update(attrs)
    
    #arr_node.flush()
    #h5f.close()
    
def write_sp(sp_dict, fname, dataset,overwrite=False):
    
    h5f = _open_file(fname, 'a')
    
    sp = sp_dict['data']
    
    parts = dataset.split('/')
    group = '/'.join(parts[:-1])
    node_name = parts[-1]

    if overwrite:
        try:
            h5f.removeNode(group, node_name)
        except tables.exceptions.NodeError:
            pass
    
    atom = tables.Atom.from_dtype(sp.dtype)
    shape = sp.shape
    filters = tables.Filters(complevel=0, complib='zlib')
    arr_node = h5f.createCArray(group, node_name, atom, shape, 
                                     filters=filters, 
                                     createparents=True)
    arr_node[:] = sp
    
    arr_node.attrs['sampfreq']=sp_dict['FS']
    #attrs = sp_dict.copy()
    #del attrs['data']
    #del attrs['FS']

    #for k, v in attrs.items():
    #    arr_node.setAttr(k, v)

    #arr_node._v_attrs.__dict__.update(attrs)
    
    #arr_node.flush()
    #h5f.close()

