#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import extract, cluster

def snr_spike(spike_waves, scale=5.):
    """Estimate signal-to-noise ratio (SNR) as a ratio of
    peak-to-peak amplitude of an average spike to the std. deviation
    of residuals"""
    
    sp_data = spike_waves['data']
    avg_spike = sp_data.mean(1)

    peak_to_peak = avg_spike.max()-avg_spike.min()

    residuals = sp_data-avg_spike[:, np.newaxis]

    noise_std = np.sqrt(residuals.var())

    snr = peak_to_peak/(noise_std*scale)

    return snr

def snr_clust(spike_waves, noise_waves):
    """Calculate signal-to-noise ratio by comparing average P2P
    amplitude of spike cluster to noise cluster.
    
    See also: extract_noise_cluster"""

    def _calc_p2p(data):
        p2p = data.max(0)-data.min(0)
        return p2p.mean()

    sp_data = spike_waves['data']
    avg_p2p_spk = _calc_p2p(sp_data)

    noise_data = noise_waves['data']
    avg_p2p_ns = _calc_p2p(noise_data)

    snr = avg_p2p_spk/avg_p2p_ns

    return snr

def extract_noise_cluster(sp, spt, sp_win, type="positive"):

    spike_waves = extract.extract_spikes(sp, spt, sp_win)
    
    if type == "positive":
        threshold = calc_noise_threshold(spike_waves, 1)
        spt_noise = extract.detect_spikes(sp, threshold, 'rising')
        spt_noise = extract.align_spikes(sp, spt_noise, sp_win, 'max')
    else:
        threshold = calc_noise_threshold(spike_waves, -1)
        spt_noise = extract.detect_spikes(sp, threshold, 'falling')
        spt_noise = extract.align_spikes(sp, spt_noise, sp_win, 'min')

    spt_noise = extract.remove_spikes(spt_noise, spt, sp_win)
    sp_waves = extract.extract_spikes(sp, spt_noise, sp_win)

    return sp_waves

def calc_noise_threshold(spike_waves, sign=1, frac_spikes=0.02, frac_max=0.5):
    """ Find threshold to extract noise cluster according to algorithm described in
    Joshua et al. (2007)

    Arguments:

    * spike_waves : dictionary
      waveshapes of spikes from the identified single unit (extracted
      with extract.extrac_spikes)

    * sign : int
      sign should be negative for negative-going spikes and positive for
      postitive-going spikes. Note that only sign of this number is taken into
      account.

    * frac_spikes : float, default 0.02
      fraction of largest (smallest) spikes to calculate the threshold
      from

    * frac_max : float
      fraction of the average peak amplitude to use as a treshold


    Returns: float

      threshold to obtain a noise cluster 

    """

    gain = np.sign(sign)

    peak_amp = np.max(gain*spike_waves['data'],0)     
    frac_spikes = 0.02
    frac_max = 0.5
    peak_amp.sort()
    threshold = frac_max*np.mean(peak_amp[:int(frac_spikes*len(peak_amp))])
    threshold *= gain

    return threshold

def calc_isolation_score(sp, spt, sp_win, spike_type='positive',
        lam=10.):
    """Calculate isolation index according to Joshua et al. (2007)
    
    Arguments:
    
    * sp : dict
      raw extracellular recordings
      
    * spt : dict
      spike times
    
    * sp_win : list or tuple
      window used for spike extraction
      
    * spike_type : positive or negative
      indicates if the spikes occuring at time points given by spt are
      positive or negative going

    * lambda : float, (0, Inf)
      determines the "softness" of clusters; 
      
    Returns:
        
      isolation_score : float
      a value from the range [0,1] indicating the quality of sorting
      (1=ideal isolation of spikes)
    """
   
    #Extract waveshapes of isolated spikes and noise
    spike_waves = extract.extract_spikes(sp, spt, sp_win)
    noise_waves = extract_noise_cluster(sp, spt, sp_win, spike_type) 
    n_spikes = spike_waves['data'].shape[1]
    
    #calculate distance between spikes and all other events
    all_waves = {'data': np.hstack((spike_waves['data'],
                                    noise_waves['data']))}
    dist_matrix = cluster.dist_euclidean(spike_waves, all_waves)
    d_0 = dist_matrix[:,:n_spikes].mean()
   
    #normalize distance and convert it to similarity
    similarity = np.exp(-dist_matrix*lam/d_0)
    similarity[dist_matrix==0]=0.
    
    #normalize probablity to 1 (softmax)
    similarity = similarity/similarity.sum(1)[:, np.newaxis]
    
    p_x = similarity[:,:n_spikes].sum(1)

    isolation_score = p_x.mean()

    return isolation_score

