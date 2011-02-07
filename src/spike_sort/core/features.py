#!/usr/bin/env python
#coding=utf-8
"""
Provides functions to calculate spike waveforms features.

Functions starting with `fet` implement various features calculated
from the spike waveshapes. They have usually one required argument
:ref:`spikewave` structure (but there are exceptions!).  

Each of the function returns a (mapping) object with following keys:

 * data -- an array of shape (n_spikes, n_features)
 * names -- a list of length n_features with feature labels
"""


import numpy as np
import matplotlib.pyplot as plt

def split_cells(features, idx, which='all'):
    """return the spike features splitted into separate cells"""

    if which == 'all':
        classes = np.unique(idx)
    else:
        classes = which
    
    data = features['data']
    names = features['names']
    feature_dict = dict([(cl, {'data': data[idx==cl, :],
                               'names': names}) for cl in classes])

    return feature_dict


def select(features_dict, features_ids):

    """Choose selected features from the collection"""
    
    def _feat2idx(id):
        if type(id) is str:
            i = np.nonzero(names==id)[0][0]
            return i
        else:
            return id
   
    names = features_dict['names']
    features = features_dict['data']
    ii = np.array([_feat2idx(id) for id in features_ids])

    selected = {"data": features[:,ii], "names": names[ii]}

    return selected

def combine(args, normalize=True):
    """Combine features into a single structure
    
    :arguments:
        * args -- a tuple of features
    
    :output:
        * combined_fetures -- dictionary with two keys `data` -- feature
          array of shape `(n_spikes, n_features)` and `names` -- feature labels
    """

    features = [d['data'] for d in args]
    names = [d['names'] for d in args]

    data  = np.hstack(features)
    if normalize:
        data = (data-data.min(0)[np.newaxis,:])
        data = data/data.max(0)[np.newaxis, :]

    combined_features = {"data": data,
                "names":np.concatenate(names)}
    return combined_features

def PCA(data,ncomps=2):
    """
    Perfrom a principle component analysis on `data` and project
    data on `ncomps` eigenvectors

    :arguments:
     
     * data -- (n_vars, n_obs) array where `n_vars` is the number of
       variables (vector dimensions) and `n_obs` the number of
       observations

    :output:

     * evals -- sorted eigenvalues
     * evecs -- sorted eigenvectors
     * score -- projection of the data on `ncomps` components
     """

    #norm=data/np.std(data,1)[:,np.newaxis]
    #norm[np.isnan(norm)]=0
    #norm = data
    data = data.astype(np.float64)
    K=np.cov(data)
    evals,evecs=np.linalg.eig(K)
    order=np.argsort(evals)[::-1]
    evecs=np.real(evecs[:,order])
    evals=np.abs(evals[order])
    score= np.dot(evecs[:,:ncomps].T,data)
    score = score/np.sqrt(evals[:ncomps, np.newaxis])
    return evals,evecs,score

def _get_data(spk_dict, contacts):
    spikes = spk_dict["data"]
    if not contacts=="all":
        contacts = np.asarray(contacts)
        try:
            spikes = spikes[:,:,contacts]
        except IndexError:
            raise IndexError("contacts must be either 'all' or a list of valid"+
                             " contact indices" )
    return spikes

def fetPCs(spikes_data,ncomps=2, contacts='all'):
    """Calculate principal components (PCs).
    
    :arguments:
        
     * spikes -- spikewave structures
     * ncomps -- number of components to retain
     
     :output:

     * pcs -- projection scores of size `(n_contacts*ncomps, n_spikes)`
     * names -- feature labels ("Ch0:PC0', "Ch0:PC1", "Ch1:PC0", etc.)
    """

    spikes = _get_data(spikes_data, contacts)
    def _getPCs(data):
        _, _, sc=PCA(data[:,:],ncomps)
        return sc

    if spikes.ndim==3:
        n_channels = spikes.shape[2]
        sc=[]
        for i in range(n_channels):
            sc+=[_getPCs(spikes[:,:,i])]
        sc=np.vstack(sc)
    else:
        sc=_getPCs(spikes)
        n_channels = 1
    sc=np.asarray(sc).T
    
    names = ["Ch%d:PC%d" % (j,i) for i in range(ncomps) for j in
            range(n_channels)]
    
    return {'data': sc, "names":names}

def fetP2P(spikes_data, contacts='all'):
    """Calculate peak-to-peak amplitudes of spike waveforms.

    :arguments:
     
     * spikes -- spikewave structure with spike waveshapes (see
       documentation for detailed specification)

    :output:

     * p2p (int) -- 2D array of peak-to-peak amplitudes in subsequent
       channels (contacts)
     
     * name -- feature labels (ex. Ch0:P2P, Ch1:P2P)

    **Example**

     We will generate a spikewave structure containing only a single
     spike on a single channel

     >>> import numpy as np
     >>> from spike_sort import features
     >>> time = np.arange(0,2*np.pi,0.01) 
     >>> spikes = np.sin(time)[:,np.newaxis, np.newaxis]
     >>> spikewave = {"data": spikes, "time":time, "contacts":1, "FS":1}
     >>> p2p = features.fetP2P(spikewave)
     >>> print p2p['data']
     [[ 1.99999683]]

    """

    spikes = _get_data(spikes_data, contacts)
            
    p2p=spikes.max(axis=0)-spikes.min(axis=0)
    if p2p.ndim<2:
        p2p=p2p[:,np.newaxis]

    names = ["Ch%d:P2P" % i for i in range(p2p.shape[1])]

    return {'data':p2p, 'names':names}

def fetSpIdx(spikes_data):
    """
    Spike sequential index (0,1,2, ...)
    """

    spikes = _get_data(spikes_data, [0])

    n_datapts = spikes.shape[1]

    return {'data':np.arange(n_datapts)[:, np.newaxis],'names':["SpIdx"]}

def fetSpTime(spt_dict):
    """
    Spike occurance time im milliseconds.

    :arguments:
        * spt_dict -- dictionary with `data` key containing spike times
    :output:
        * spt -- spike times
        * label
    """

    spt = spt_dict['data']

    return {'data': spt[:, np.newaxis], 'names':["SpTime"]}
