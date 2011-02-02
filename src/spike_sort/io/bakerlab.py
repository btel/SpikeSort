#!/usr/bin/env python
#coding=utf-8

import os
import numpy as np
import json
import re
from tempfile import mkdtemp
import tables

_regexp="^/(?P<subject>[a-zA-z]+)/s(?P<ses_id>.+)/el(?P<el_id>[0-9]+)$"

def read_sp(conf_file, dataset, memmap=None):
    """Reads raw spike waveform from file in bakerlab format
    
    :arguments:
    - conf_file -- configuration file
    - dataset -- dataset path
    - memmap -- use memmory mapped arrays to save some memory
    """
    with file(conf_file) as fid:
        conf_dict = json.load(fid)
    m = re.match(_regexp, dataset)
    rec_dict = m.groupdict()
    n_contacts = conf_dict['n_contacts']
    f_spike = conf_dict['fspike']
    
    dirname = conf_dict['dirname']
    sp_list = []
    rec_dict['contact_id']=1
    
    full_path = os.path.join(dirname, f_spike)
    fname = full_path.format(**rec_dict)
    sp = np.fromfile(fname, dtype=np.int16)/200.

    
    if memmap=="numpy":
        #create temporary memmory mapped array
        filename = os.path.join(mkdtemp(), 'newfile.dat')
        fp = np.memmap(filename, dtype='float', mode='w+', 
                       shape=(len(sp),n_contacts))
    elif memmap=="tables":
        atom = tables.Atom.from_dtype(sp.dtype)
        shape = (len(sp), n_contacts)
        filters = tables.Filters(complevel=3, complib='blosc')
        filename = os.path.join(mkdtemp(), 'newfile.dat')
        h5f = tables.openFile(filename,'w')
        fp = h5f.createCArray('/', "test", atom, shape, filters=filters)
    else:
        fp = np.empty((len(sp), n_contacts), dtype='float')

    fp[:,0]=sp
    for i in range(1,n_contacts):
        rec_dict['contact_id']=i+1
        fname = full_path.format(**rec_dict)
        sp = np.fromfile(fname,dtype=np.int16)
        fp[:,i]=sp/200.
    del sp
    return {'data':fp, "FS":conf_dict['FS'], "n_contacts":n_contacts} 

def write_sp(sp_dict, conf_file, dataset):
    """Write raw spike waveform to a file in bakerlab format
    
    :arguments:
    - sp_dict -- spike waveform dict
    - conf_file -- configuration file
    - dataset -- dataset path
    """
    sp = sp_dict['data']
    with file(conf_file) as fid:
        conf_dict = json.load(fid)
    m = re.match(_regexp, dataset)
    n_contacts = sp.shape[1]
    rec_dict = m.groupdict()
    
    for i in range(n_contacts):
        rec_dict['contact_id']=i+1
        fname = conf_dict['fspike'].format(**rec_dict)
        full_path = os.path.join(conf_dict['dirname'], fname)
        sp_int = (sp[:,i]*200).astype(np.int16)
        sp_int.tofile(full_path)

def read_spt(dir_name, dataset):
    """Returns spike times in miliseconds:
    
    :arguments:
    -- dir_name : directory names with the data
    -- dataset : dataset path
    """
    
    fname = os.path.join(dir_name, dataset+".spt")
    spt = np.fromfile(fname, dtype=np.int32)
    return {"data": spt/200.}


def write_spt(spt_dict, dir_name, dataset):
    """Returns spike times in miliseconds:
    
    Arguments:
    -- dir_name : directory names with the data
    -- dataset : dataset name
    """

    spt = spt_dict['data']
    
    fname = os.path.join(dir_name, dataset+".spt")
    export_spt = (spt*200).astype(np.int32)
    export_spt.tofile(fname)
