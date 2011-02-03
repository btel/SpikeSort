#!/usr/bin/env python
#coding=utf-8

import numpy as np
from scipy import interpolate

def split_cells(spikes, idx, which='all'):
    """return the spike features splitted into separate cells"""

    if which == 'all':
        classes = np.unique(idx)
    else:
        classes = which
    
    data = spikes['data']
    time = spikes['time']
    spikes_dict = dict([(cl, {'data': data[:,idx==cl, :], 'time': time}) 
                        for cl in classes])

    return spikes_dict

def remove_spikes(spt_dict, remove_dict, tolerance):
    spt_data = spt_dict['data']
    spt_remove = remove_dict['data']

    min, max = tolerance

    for t in spt_remove:
        spt_data = spt_data[(spt_data>(t+max)) | (spt_data<(t+min))]

    spt_ret = spt_dict.copy()

    spt_ret['data'] = spt_data

    return spt_ret

def detect_spikes(spike_data, thresh='auto', edge="rising",
                  contact=0):
    """Detects spikes in extracellular data using amplitude thresholding.

    Arguments:

    -- spike_data : dict
       extracellular waveforms

    -- thresh : float or 'auto'
       threshold for detection. if thresh is 'auto' it will be
       estimated from the data.

    -- edge : {"rising" or "falling"}
       which edge to trigger on

    -- contact: index of tetrode contact to use for detection

    Returns: 
    dictionary with 'data' key which contains detected threshold
    crossing in miliseconds

    """
    
    sp_data = spike_data['data']
    n_contacts = spike_data['n_contacts']

    #if n_contacts>1:
    #    sp_data = sp_data[:,contact]

    FS = spike_data['FS']

    if thresh=='auto':
        thresh = 8*np.sqrt(float(np.var(sp_data[:10*FS,contact])))
        if edge == 'falling':
            thresh = -thresh

    if edge == "rising":
        i, = np.where((sp_data[:-1,contact]<thresh) & (sp_data[1:,contact]>thresh))
    elif edge == "falling":
        i, = np.where((sp_data[:-1,contact]>thresh) & (sp_data[1:,contact]<thresh))
    else:
        raise "Edge must be 'rising' or 'falling'"
    spt = i*1000./FS

    spt_dict = {'data': spt}

    return spt_dict

def extract_spikes(spike_data, spt_dict, sp_win, resample=None,
                   contacts='all'):
    """Returns spike wave shapes.

    Arguments:

    -- spike_data: extracelluler waveforms (n_pts, n_spikes, n_contacts)
    -- spt : spike times
    -- sp_win : temporal extent of the wave shape 

    """

    sp_data = spike_data['data']
    n_contacts = spike_data['n_contacts']

    if contacts == "all":
        contacts = np.arange(n_contacts)
    elif type(contacts) is int:
        contacts = np.array([contacts])
    else:
        contacts = np.asarray(contacts)

    FS = spike_data['FS']
    spt = spt_dict['data']

    indices = (spt/1000.*FS).astype(np.int32)
    win = (np.asarray(sp_win)/1000.*FS).astype(np.int32)
   
    isOutOfBound = ((indices>(sp_data.shape[0]-win[1])) | (indices<-win[0]))
    correct_indices = indices[~isOutOfBound]
    truncated_indices = indices[isOutOfBound]

    time = np.arange(win[1]-win[0])*1000./FS+sp_win[0]


    if resample is None or resample==1:
        spWave = np.zeros((len(time), len(spt), len(contacts)), dtype=np.float32)
        for i,sp in enumerate(correct_indices):
            spWave[:,i,:] = sp_data[sp+win[0]:sp+win[1],contacts]
        return {"data":spWave, "time": time, "FS": FS}

    else:
        FS_new = FS*resample
        resamp_time = np.arange(sp_win[0], sp_win[1], 1000./FS_new)
        spWave = np.zeros((len(resamp_time), len(spt), len(contacts)), dtype=np.float32)
       
        for i,sp in enumerate(correct_indices):
            time = np.arange(sp+win[0]-1, sp+win[1]+1)*1000./FS
            for contact in contacts:
                new_wave  = sp_data[sp+win[0]-1:sp+win[1]+1, contact]
                tck = interpolate.splrep(time, new_wave, s=0)
                spWave[:,i,contact] = interpolate.splev(resamp_time+spt[i], tck, der=0)

        return {"data":spWave, "time": resamp_time, "FS": FS_new}

def resample_spikes(spikes_dict, FS_new):

    sp_waves = spikes_dict['data']
    time = spikes_dict['time']
    FS = spikes_dict['FS']

    resamp_time = np.arange(time[0], time[-1], 1000./FS_new)
    n_pts, n_spikes, n_contacts = sp_waves.shape

    spike_resamp = np.empty((len(resamp_time), n_spikes, n_contacts))

    for i in range(n_spikes):
        for contact in range(n_contacts):
            tck = interpolate.splrep(time, sp_waves[:, i, contact],s=0)
            spike_resamp[:,i, contact] = interpolate.splev(resamp_time, tck, der=0)

    return {"data":spike_resamp, "time":resamp_time, "FS":FS}
    


def align_spikes(spike_data, spt_dict, sp_win, type="max", resample=1,
                contact=0):
    
    """aligns spike waves and returns corrected spike times"""

    spt = spt_dict['data']
    
    sp_waves_dict = extract_spikes(spike_data, spt_dict, sp_win,
            resample=resample, contacts=contact)

    sp_waves = sp_waves_dict['data'][:,:,0]
    time = sp_waves_dict['time']

    if type=="max":
        i = sp_waves.argmax(0)
        spt_new = spt + time[i]
    elif type=="min":
        i = sp_waves.argmin(0)
        spt_new = spt + time[i]


    return {"data": spt_new}


def merge_spikes(spike_waves1, spike_waves2):
    """Merges two sets of spike waves
    
    Arguments:
    
    * spike_waves1 : dictionary
    * spike_waves2 : dictionary
      both spike wave sets must be defined within the same time window
      and with the same sampling frequency

    Returns:

    * spike_waves : dictionary
      merged spike waveshapes

    * clust_idx : array
      labels denoting to which set the given spike originally belonged
      to
    
    """

    sp_data1 = spike_waves1['data']
    sp_data2 = spike_waves2['data']

    sp_data = np.hstack((sp_data1, sp_data2))
    spike_waves = spike_waves1.copy()
    spike_waves['data'] = sp_data

    clust_idx = np.concatenate((np.ones(sp_data1.shape[1]),
                               np.zeros(sp_data2.shape[1])))

    return spike_waves, clust_idx

def merge_spiketimes(spt1, spt2, sort=True):
    """Merges two sets of spike times 
    
    Arguments:
    
    * spt1 : dictionary
    * spt2 : dictionary
    * sort : bool
      False if you don't want to be the spike times sorted.

    Returns:

    * spt : dictionary
      dictionary with merged spike time arrrays under data key

    * clust_idx : array
      labels denoting to which set the given spike originally belonged
      to
    
    """

    spt_data1 = spt1['data']
    spt_data2 = spt2['data']

    spt_data = np.concatenate((spt_data1, spt_data2))
    i = spt_data.argsort()
    spt_data = spt_data[i]


    clust_idx = np.concatenate((np.ones(spt_data1.shape[0]),
                               np.zeros(spt_data2.shape[0])))
    clust_idx = clust_idx[i]
    spt_dict = {"data": spt_data}

    return spt_dict, clust_idx
