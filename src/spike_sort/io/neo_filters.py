#!/usr/bin/env python
#coding=utf-8

import spike_sort as sort
from spike_beans import components
from spike_beans.components import GenericSource
try:
    import neo
except ImportError:
    raise ImportError, "To use the extra data filters you have to install Neo package"
import numpy as np
import os

class AxonFilter(object):
    """Read Axon .abf files
    
    Parameters
    ----------
    fname : str
        path to file

    electrodes : (optional) list of ints
        list of electrode indices to use
    
    """
   
    def __init__(self, fname, electrodes=None):
        self.reader = neo.io.AxonIO(fname)
        self.block = self.reader.read_block()
        self.electrodes = electrodes

    def read_sp(self, dataset=None):
        electrodes = self.electrodes
        analogsignals = self.block.segments[0].analogsignals
        if electrodes is not None:
            analogsignals = [analogsignals[i] for i in electrodes]
        sp_raw = np.array(analogsignals)
        FS = float(analogsignals[0].sampling_rate.magnitude)
        n_contacts, _ = sp_raw.shape
        return {"data": sp_raw, "FS": FS, "n_contacts": n_contacts}
    
    def write_sp(self):
        raise NotImplementedError, "Writing to Axon files not yet implemented"
    def write_spt(self, spt_dict, dataset, overwrite=False):
        raise NotImplementedError, "Writing spikes in Axon format not yet implemented"

class NeoSource(components.GenericSource):
    def __init__(self, fname, electrodes=None, overwrite=False):
        GenericSource.__init__(self, '', overwrite)
        
        root, ext = os.path.splitext(fname)
        ext = ext.lower()
        if ext == '.abf':
            self.io_filter = AxonFilter(fname, electrodes)
        else:
            raise IOError, "Format {0} not recognised".format(ext)
        self.read_sp = self.io_filter.read_sp
        self.write_sp = self.io_filter.write_sp
        self.write_spt = self.io_filter.write_spt
