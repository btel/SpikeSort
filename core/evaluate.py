#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
import extract, cluster
import warnings

def deprecation(message):
    warnings.warn(message, DeprecationWarning, stacklevel=2)

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


    deprecation("extract_noise_cluster deprecated. Use"
                " detect_noise and extract_spikes instead.")
    spt_noise = detect_noise_spikes(sp, spt, sp_win, type)
    sp_waves = extract.extract_spikes(sp, spt_noise, sp_win)

    return sp_waves

def detect_noise(sp, spt, sp_win, type="positive"):

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

    return spt_noise

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

def isolation_score(sp, spt, sp_win, spike_type='positive', lam=10., max_spikes=None):
    
    spike_waves = extract.extract_spikes(sp, spt, sp_win)
    spt_noise = detect_noise(sp, spt, sp_win, spike_type)
    noise_waves = extract.extract_spikes(sp, spt_noise, sp_win)

    iso_score = calc_isolation_score(spike_waves, noise_waves,
            spike_type, lam=lam, max_spikes=max_spikes)

    return iso_score

def _iso_score_dist(dist, lam, n_spikes):

    """Calculate isolation score from a distance matrix

    Arguments:
    * dist : numpy array
      NxM matrix, where N is number of spikes and M is number of all
      events (spikes + noise)

    * lam : float 
      lambda parameter

    * n_spikes : int
      number of spikes (N)

    Returns:
    * isolation score
    """


    distSS = dist[:, :n_spikes]
    distSN = dist[:, n_spikes:] 

    d0 = distSS.mean()
    expSS = np.exp(-distSS*lam*1./d0);
    expSN = np.exp(-distSN*lam*1./d0);
 
    sumSS = np.sum(expSS - np.eye(n_spikes),1)
    sumSN = np.sum(expSN, 1) 

    correctProbS = sumSS /  (sumSS + sumSN)
    isolation_score = correctProbS.mean()

    return isolation_score



def calc_isolation_score(spike_waves, noise_waves, spike_type='positive',
        lam=10., max_spikes=None):
    """Calculate isolation index according to Joshua et al. (2007)
    
    Arguments:
    
    * spike_waves : dict
      
    * noise_waves : dict
    
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

    #Memory issue: sample spikes if too many
    if max_spikes is not None:
        if spike_waves['data'].shape[1]>max_spikes:
           i = np.random.rand(max_spikes).argsort()
           spike_waves = spike_waves.copy()
           spike_waves['data'] = spike_waves['data'][:, i]
        if noise_waves['data'].shape[1]>max_spikes:
           i = np.random.rand(max_spikes).argsort()
           noise_waves = noise_waves.copy()
           noise_waves['data'] = noise_waves['data'][:, i]


    n_spikes = spike_waves['data'].shape[1]
    
    #calculate distance between spikes and all other events
    all_waves = {'data': np.hstack((spike_waves['data'],
                                    noise_waves['data']))}
    dist_matrix = cluster.dist_euclidean(spike_waves, all_waves)
    #d_0 = dist_matrix[:,:n_spikes].mean()
   
    isolation_score = _iso_score_dist(dist_matrix, lam,
            n_spikes)

    return isolation_score

