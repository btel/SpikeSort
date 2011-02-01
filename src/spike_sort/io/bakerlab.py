#!/usr/bin/env python
#coding=utf-8

import os
import numpy as np
import json
import re

_regexp="^/(?P<subject>[a-zA-z]+)/s(?P<ses_id>.+)/el(?P<el_id>[0-9]+)$"
def read_sp(conf_file, dataset):
    with file(conf_file) as fid:
        conf_dict = json.load(fid)
    m = re.match(_regexp, dataset)
    fname = conf_dict['fspike'].format(**m.groupdict())
    full_path = os.path.join(conf_dict['dirname'], fname)
    sp = np.fromfile(full_path,dtype=np.int16)/200.
    return {'data':sp, "FS":conf_dict['FS']} 

def write_sp(sp_dict, conf_file, dataset):
    sp = sp_dict['data']
    with file(conf_file) as fid:
        conf_dict = json.load(fid)
    m = re.match(_regexp, dataset)
    fname = conf_dict['fspike'].format(**m.groupdict())
    full_path = os.path.join(conf_dict['dirname'], fname)
    sp_int = (sp*200).astype(np.int16)
    sp_int.tofile(full_path)

def read_spt(dir_name, dataset):
    """Returns spike times in miliseconds:
    
    Arguments:
    -- dir_name : directory names with the data
    -- dataset : dataset name
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
