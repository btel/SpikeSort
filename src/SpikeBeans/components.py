'''
Created on Feb 9, 2011

@author: bartosz
'''

import spike_sort as sort

import base
from spike_sort.io.filters import BakerlabFilter, PyTablesFilter
from spike_sort import features
import numpy as np

class BakerlabSource(base.Component, BakerlabFilter):  
    
    def __init__(self, conf_file, dataset):
        BakerlabFilter.__init__(self, conf_file)
        self.dataset = dataset
        self._signal = None
        self._events = None
        super(BakerlabSource, self).__init__()
        
    def read_signal(self):
        if self._signal is None:
            self._signal = self.read_sp(self.dataset)
        return self._signal
    
    def read_events(self):
        if self._events is None:
            self._events = self.read_spt(self.dataset)
        return self._events
    
    signal = property(read_signal)
    events = property(read_events)
    

class SpikeDetector(base.Component):
    """Detect Spikes with alignment"""
    waveform_src = base.RequiredFeature("SignalSource", 
                                        base.HasAttributes("signal"))

    def __init__(self, thresh='auto', 
                 contact=0, 
                 type='max', 
                 resample=1, 
                 sp_win=(-0.2, 0.8)):
        self.thresh = thresh
        self.contact = contact
        self.type = "max"
        self.resample = resample
        self.sp_win = sp_win
        self.sp_times = None
        super(SpikeDetector, self).__init__()
        
    def _detect(self):
        sp = self.waveform_src.signal
        spt = sort.extract.detect_spikes(sp,  contact=self.contact,
                                               thresh=self.thresh)
        self.sp_times = sort.extract.align_spikes(sp, spt, 
                                                  self.sp_win, 
                                                  type=self.type, 
                                                  resample=self.resample)
    
    def update(self):
        self._detect()
        super(SpikeDetector, self).update()

    def read_events(self):
        if self.sp_times is None:
            self._detect()
        return self.sp_times
    
    events = property(read_events)
    
class SpikeExtractor(base.Component):
    waveform_src = base.RequiredFeature("SignalSource", 
                                    base.HasAttributes("signal"))
    spike_times = base.RequiredFeature("SpikeMarkerSource", 
                                    base.HasAttributes("events"))
    
    def __init__(self, sp_win=[-0.2,0.8]):
        self._sp_shapes = None
        self.sp_win = sp_win
        super(SpikeExtractor, self).__init__()
    
    def _extract_spikes(self):
        sp = self.waveform_src.signal
        spt = self.spike_times.events
        self._sp_shapes = sort.extract.extract_spikes(sp, spt, self.sp_win)
    
    def read_spikes(self):
        if self._sp_shapes is None:
            self._extract_spikes()
        return self._sp_shapes
    
    def update(self):
        self._extract_spikes()
        super(SpikeExtractor, self).update()
    
    spikes = property(read_spikes)
    
class FeatureExtractor(base.Component):
    spikes_src = base.RequiredFeature("SpikeSource", 
                                      base.HasAttributes("spikes"))
    
    def __init__(self, normalize=True):
        self.feature_methods = []
        self._feature_data = None
        self.normalize = normalize
        super(FeatureExtractor, self).__init__()
        
    def add_feature(self, name, *args, **kwargs):
        func_name = "fet" + name
        func = lambda x: features.__getattribute__(func_name)(x, *args, **kwargs)
        self.feature_methods.append(func)
    
    def _calc_features(self):
        spikes = self.spikes_src.spikes
        feats = [f(spikes) for f in self.feature_methods]
        self._feature_data = features.combine(feats, normalize=self.normalize)
    
    def read_features(self):
        if self._feature_data is None:
            self._calc_features()
        return self._feature_data

    def update(self):
        self._calc_features()
        super(FeatureExtractor, self).update()
    
    features = property(read_features)
        
class ClusterAnalyzer(base.Component):
    feature_src = base.RequiredFeature("FeatureSource", 
                                       base.HasAttributes("features"))
    
    def __init__(self, method, *args, **kwargs):
        self.method = method
        self.args = args
        self.kwargs = kwargs
        self.cluster_labels = None
        super(ClusterAnalyzer, self).__init__()
    
    def _cluster(self, idx, method, *args,**kwargs):
        feature_data = self.feature_src.features
        if idx is not None:
            feature_data['data'] = feature_data['data'][idx,:]
            clust_idx = sort.cluster.cluster(method, feature_data, *args, 
                                         **kwargs)
            all_labels = set(range(1, 100))
            used_labels = set(np.unique(self.cluster_labels))
            free_labels = list(all_labels - used_labels)
            free_labels.sort(reverse=True)
            for l in np.unique(clust_idx):
                new_label = free_labels.pop()
                clust_idx[clust_idx==l]= new_label
                
            self.cluster_labels[idx] = clust_idx
        else:
            clust_idx = sort.cluster.cluster(method,feature_data, *self.args, 
                                         **kwargs)
            self.cluster_labels = clust_idx+1
            
    
    def read_labels(self):
        if self.cluster_labels is None:
            self._cluster(None, self.method, *self.args, **self.kwargs)
        return self.cluster_labels
    
    def recluster(self, label, method=None, *args, **kwargs):
        if method is None:
            method = self.method
        if not args:
            args = self.args
        if not kwargs:
            kwargs = self.kwargs 
        self._cluster(self.cluster_labels==label, method, *args, **kwargs)
        
    def delete_cell(self, cell_id):
        self.cluster_labels[self.cluster_labels==cell_id] = 0
    
    def merge_cells(self, cell_ids):
        for cell in cell_ids:
            self.cluster_labels[self.cluster_labels==cell]=cell_ids[0]
        
    def update(self):
        self._cluster(None, self.method, *self.args, **self.kwargs)
        super(ClusterAnalyzer, self).update()
            
    labels = property(read_labels)
                
        
        
            