#!/usr/bin/env python
#coding=utf-8

import os
import numpy as np

def read_spt(dir_name, dataset):
    """Returns spike times in miliseconds:
    
    Arguments:
    -- dir_name : directory names with the data
    -- dataset : dataset name
    """
    
    fname = os.path.join(dir_name, dataset+".spt")
    spt = np.fromfile(fname, dtype=np.int32)
    return spt/200.


def write_spt(spt, dir_name, dataset):
    """Returns spike times in miliseconds:
    
    Arguments:
    -- dir_name : directory names with the data
    -- dataset : dataset name
    """
    
    fname = os.path.join(dir_name, dataset+".spt")
    export_spt = (spt*200).astype(np.int32)
    export_spt.tofile(fname)
