#!/usr/bin/env python
#coding=utf-8

import os
import numpy as np
import json
import re

_regexp="^/(?P<subject>[a-zA-z]+)/s(?P<ses_id>.+)/el(?P<el_id>[0-9]+)$"

def read_sp(conf_file, dataset):
    """Reads raw spike waveform from file in bakerlab format
    
    :arguments:
    - conf_file -- configuration file
    - dataset -- dataset path
    """
    with file(conf_file) as fid:
        conf_dict = json.load(fid)
    m = re.match(_regexp, dataset)
    rec_dict = m.groupdict()
    n_contacts = conf_dict['n_contacts']
    f_spike = conf_dict['fspike']
    dirname = conf_dict['dirname']
    sp_list = []
    for i in range(n_contacts):
        rec_dict['contact_id']=i+1
        fname = f_spike.format(**rec_dict)
        full_path = os.path.join(dirname, fname)
        sp = np.fromfile(full_path,dtype=np.int16)/200.
        sp_list.append(sp)
    sp_list = np.array(sp_list).T
    return {'data':sp_list, "FS":conf_dict['FS']} 

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
