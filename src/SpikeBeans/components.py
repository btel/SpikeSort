'''
Created on Feb 9, 2011

@author: bartosz
'''

import spike_sort as sort

import base
from spike_sort.io.filters import BakerlabFilter, PyTablesFilter
from spike_sort import features

class BakerlabSource(base.Component, BakerlabFilter):
    
    def __init__(self, conf_file, dataset):
        BakerlabFilter.__init__(self, conf_file)
        self.dataset = dataset
        self.signal = None
        self.events = None
        
    def read_signal(self):
        if self.signal is None:
            self.signal = self.read_sp(self.dataset)
        return self.signal
    
    def read_events(self):
        if self.events is None:
            self.events = self.read_spt(self.dataset)
        return self.events
    

class SpikeDetector(base.Component):
    """Detect Spikes with alignment"""
    waveform_src = base.RequiredFeature("SignalSource", 
                                        base.HasAttributes("read_signal"))

    def __init__(self, thresh='auto', 
                 contact=0, 
                 type='max', 
                 resample=1, 
                 sp_win=[-0.2, 0.8]):
        self.thresh = thresh
        self.contact = contact
        self.type = "max"
        self.resample = resample
        self.sp_win = sp_win
        self.sp_times = None
        
    def _detect(self):
        sp = self.waveform_src.read_signal()
        spt = sort.extract.detect_spikes(sp,  contact=self.contact,
                                               thresh=self.thresh)
        self.sp_times = sort.extract.align_spikes(sp, spt, 
                                                  self.sp_win, 
                                                  type=self.type, 
                                                  resample=self.resample)
        
    def read_events(self):
        if self.sp_times is None:
            self._detect()
        return self.sp_times
    
class SpikeExtractor(base.Component):
    waveform_src = base.RequiredFeature("SignalSource", 
                                    base.HasAttributes("read_signal"))
    spike_times = base.RequiredFeature("SpikeMarkerSource", 
                                    base.HasAttributes("read_events"))
    
    def __init__(self, sp_win=[-0.2,0.8]):
        self.sp_shapes = None
        self.sp_win = sp_win
    
    def _extract_spikes(self):
        sp = self.waveform_src.read_signal()
        spt = self.spike_times.read_events()
        self.sp_shapes = sort.extract.extract_spikes(sp, spt, self.sp_win)
    
    def read_spikes(self):
        if self.sp_shapes is None:
            self._extract_spikes()
        return self.sp_shapes
    
class FeatureExtractor(base.Component):
    spikes_src = base.RequiredFeature("SpikeSource", 
                                      base.HasAttributes("read_spikes"))
    
    def __init__(self, normalize=True):
        self.features = []
        self.feature_data = None
        self.normalize = normalize
        
    def add_feature(self, name, *args, **kwargs):
        func_name = "fet" + name
        func = lambda x: features.__getattribute__(func_name)(x, *args, **kwargs)
        self.features.append(func)
    
    def _calc_features(self):
        spikes = self.spikes_src.read_spikes()
        feats = [f(spikes) for f in self.features]
        self.feature_data = features.combine(feats, normalize=self.normalize)
    
    def read_features(self):
        if self.feature_data is None:
            self._calc_features()
        return self.feature_data
        
class ClusterAnalyzer(base.Component):
    feature_src = base.RequiredFeature("FeatureSource", 
                                       base.HasAttributes("read_features"))
    
    def __init__(self, method, *args, **kwargs):
        self.method = method
        self.args = args
        self.kwargs = kwargs
        self.cluster_labels = None
    
    def _cluster(self, idx, method, *args,**kwargs):
        feature_data = self.feature_src.read_features()
        if idx is not None:
            feature_data['data'] = feature_data['data'][idx,:]
            clust_idx = sort.cluster.cluster(method, feature_data, *args, 
                                         **kwargs)
            self.cluster_labels[idx] = clust_idx
        else:
            clust_idx = sort.cluster.cluster(method,feature_data, *self.args, 
                                         **kwargs)
            self.cluster_labels = clust_idx
            
    
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
                
        
        
            