#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

def extract_spikes(spike_data, spt, sp_win):
    """Returns spike wave shapes.

    Arguments:

    -- spike_data: extracelluler waveforms
    -- spt : spike times
    -- sp_win : temporal extent of the wave shape 

    """

    sp_data = spike_data['data']
    FS = spike_data['FS']
    
    indices = (spt/1000.*FS).astype(np.int32)
    win = (np.asarray(sp_win)/1000*FS).astype(np.int32)
    noExt=win[1]-win[0]

    spWave = np.empty((noExt, len(spt)))
    time = np.arange(noExt)*1000./FS + sp_win[0] 
    for i,sp in enumerate(indices):
        spWave[:,i] = sp_data[sp+win[0]:sp+win[1]]

    return {"data":spWave, "time": time, "FS": FS}

