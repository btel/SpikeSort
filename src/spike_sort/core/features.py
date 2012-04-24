#!/usr/bin/env python
#coding=utf-8
"""
Provides functions to calculate spike waveforms features.

Functions starting with `fet` implement various features calculated
from the spike waveshapes. They have usually one required argument
:ref:`spike_wave` structure (but there are exceptions!).  

Each of the function returns a (mapping) object with following keys:

 * data -- an array of shape (n_spikes, n_features)
 * names -- a list of length n_features with feature labels
 """

import numpy as np

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

def select_spikes(features, idx):
    """Truncate features array to selected spikes. This method should be
    used to properly truncate the "features" structure
    
    Parameters
    ----------
    features : features structure
        features structure to truncate
    idx : bool list or int list
        indices of selected spikes
        
    Returns
    -------
    new_feats : features structure
        new features structure containing only data for selected spike
        indices
    """
    
    new_feats = features.copy()
    new_feats['data'] = features['data'][idx, :]
    if features.has_key('is_valid'):
        new_feats['is_valid'] = features['is_valid'][idx]
    return new_feats

def combine(args, norm=True):
    """Combine features into a single structure
    
    Parameters
    ----------
    args : tuple or list of dict
        a tuple of feature data structures
    
    Returns
    -------
    combined_fetures : dict
    """

    features = [d['data'] for d in args]
    names = [d['names'] for d in args]
    
    #get mask, if it exist
    mask = [d['is_valid'] for d in args if 'is_valid' in d]
    
    try:
        if mask:
            #combine masks using AND
            mask = reduce(np.logical_and, mask)
        data  = np.hstack(features)
    except ValueError:
        raise ValueError, 'all features must contain the same number of spikes'
    
    combined_features = {"data": data,
                         "names":np.concatenate(names)}
    if list(mask): combined_features["is_valid"] = mask
    
    if norm:
        normalize(combined_features, copy=False)

    
    return combined_features


def add_mask(feature_function):
    """Decorator to copy mask from waveshapes to features"""
    
    def _decorated(spike_data, *args, **kwargs):
        feature_data = feature_function(spike_data, *args, **kwargs)
        if 'is_valid' in spike_data:
            feature_data['is_valid'] = spike_data['is_valid']
        return feature_data
    _decorated.__doc__ = feature_function.__doc__
    return _decorated
    
    

def normalize(features, copy=True):
    """Normalize features"""
    if copy:
        features_norm = features.copy()
    else:
        features_norm = features
        
    data = features_norm['data']
    data = (data-data.min(0)[np.newaxis,:])
    data = data/data.max(0)[np.newaxis,:]
    
    features_norm['data'] = data

    return features_norm    

def PCA(data,ncomps=2):
    """Perfrom a principle component analysis.

    Parameters
    ----------
     
    data : array 
        (n_vars, n_obs) array where `n_vars` is the number of
        variables (vector dimensions) and `n_obs` the number of
        observations

    Returns
    -------
    evals : array
        sorted eigenvalues
    evecs : array 
        sorted eigenvectors
    score : array
        projection of the data on `ncomps` components
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

@add_mask
def fetPCs(spikes_data,ncomps=2, contacts='all'):
    """Calculate principal components (PCs).
    
    Parameters
    ----------
    spikes : dict
    ncomps : int, optional
        number of components to retain
     
    Returns
    -------
    features : dict
    
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

@add_mask
def fetP2P(spikes_data, contacts='all'):
    """Calculate peak-to-peak amplitudes of spike waveforms.

    Parameters
    ----------
    spikes : dict

    Returns
    -------
    features : dict

    Examples
    --------

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
    p2p = np.ptp(spikes, axis=0)
    if p2p.ndim < 2:
        p2p = p2p[:, np.newaxis]

    names = ["Ch%d:P2P" % i for i in range(p2p.shape[1])]
    return {'data': p2p, 'names': names}

@add_mask
def fetSpIdx(spikes_data):
    """Spike sequential index (0,1,2, ...)
    """

    spikes = _get_data(spikes_data, [0])

    n_datapts = spikes.shape[1]

    return {'data':np.arange(n_datapts)[:, np.newaxis],'names':["SpIdx"]}

@add_mask
def fetSpTime(spt_dict):
    """Spike occurrence time in milliseconds.
    """

    spt = spt_dict['data']

    return {'data': spt[:, np.newaxis], 'names':["SpTime"]}

@add_mask
def fetSpProjection(spikes_data, labels, cell_id=1):
    """Projection coefficient of spikes on an averaged waveform
    
    Parameters
    ----------
    spikes_data : dict
        waveform data
    labels : array
        array of length equal to number of spikes that contains 
        cluster labels
    cell_id : int
        label of cell on which all spikes should be projected.
    
    Notes
    -----
    `labels` can be also a boolean array in which case only spikes for 
    which label is True value will be averaged to determine projection
    coefficient
    """ 
    
    spikes = spikes_data["data"]
    
    proj_matrix = np.mean(spikes[:, labels==cell_id, :],1)
    projection = (proj_matrix[:,np.newaxis,:]*spikes).sum(0)
    projection /= np.sqrt((spikes**2).sum(0)*(proj_matrix**2).sum(0))
    
    ncomps = len(projection)
    names = ["Ch%d:Proj" % i for i in range(ncomps)]
    return {'data': projection, 'names':names}
    
