#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

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

def resample_spikes(spikes_dict, FS_new):

    sp_waves = spikes_dict['data']
    time = spikes_dict['time']
    FS = spikes_dict['FS']

    resamp_time = np.arange(time[0], time[-1], 1000./FS_new)
    n_spikes = sp_waves.shape[1]

    spike_resamp = np.empty((len(resamp_time), n_spikes))

    for i in range(n_spikes):
        tck = interpolate.splrep(time, sp_waves[:, i],s=0)
        spike_resamp[:,i] = interpolate.splev(resamp_time, tck, der=0)

    return {"data":spike_resamp, "time":resamp_time, "FS":FS}
    


def align_spikes(spike_data, spt, sp_win, type="max", resample=None):
    """aligns spike waves and returns corrected spike times"""


    sp_waves_dict = extract_spikes(spike_data, spt, sp_win)

    if resample:
        sp_waves_dict = resample_spikes(sp_waves_dict,
                sp_waves_dict['FS']*resample)

    sp_waves = sp_waves_dict['data']
    time = sp_waves_dict['time']

    if type=="max":
        i = sp_waves.argmax(0)
    elif type=="min":
        i = sp_waves.argmin(0)

    spt_new = spt + time[i]

    return spt_new




