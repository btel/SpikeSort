'''
Created on Feb 9, 2011

@author: bartosz
'''

import spike_sort as sort

import base
from spike_sort.io.filters import BakerlabFilter, PyTablesFilter
from spike_sort import features
from spike_sort.ui import plotting
import numpy as np

class BakerlabSource(base.Component, BakerlabFilter):  
    
    def __init__(self, conf_file, dataset, overwrite=False):
        BakerlabFilter.__init__(self, conf_file)
        self.dataset = dataset
        self._signal = None
        self._events = None
        self.overwrite = overwrite
        super(BakerlabSource, self).__init__()
        
    def read_signal(self):
        if self._signal is None:
            self._signal = self.read_sp(self.dataset)
        return self._signal
    
    def read_events(self, cell):
        node =  '/'.join((self.dataset, cell))
        if self._events is None:
            self._events = self.read_spt(node)
        return self._events
    
    def write_signal(self, signal):
        self.write_sp(signal, self.dataset)
    
    def write_events(self, cell, spt):
        node = '/'.join((self.dataset, cell))
        self.write_spt(spt, node, overwrite=self.overwrite)
    
    signal = property(read_signal, write_signal)
    events = base.dictproperty(read_events, write_events)
    
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
        self.trash_label = 0
        super(ClusterAnalyzer, self).__init__()
    
    def _cluster(self, idx, method, *args,**kwargs):
        feature_data = self.feature_src.features
        if idx is not None:
            new_features  = feature_data.copy()
            new_features['data'] = new_features['data'][idx,:]
            clust_idx = sort.cluster.cluster(method, new_features, *args, 
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
    
    def relabel(self):
        """rename cells in sequential order"""
        labels = list(np.unique(self.labels))
        if self.trash_label in labels:
            labels.remove(self.trash_label)
        for i, l in enumerate(labels):
            self.cluster_labels[self.cluster_labels==l] = i+1
        self.notify_observers()
    
    def recluster(self, label, method=None, *args, **kwargs):
        """repeat clustering of a selected cell"""
        if method is None:
            method = self.method
        if not args:
            args = self.args
        if not kwargs:
            kwargs = self.kwargs 
        self._cluster(self.cluster_labels==label, method, *args, **kwargs)
        self.notify_observers()
        
    def delete_cells(self, *cell_ids):
        """move selected labels to thrash (cluster 0). 
        if 'all' thrash all cells """
        if len(cell_ids)==1 and cell_ids[0]=='all':
            cell_id = np.unique(self.labels)
        for cell_id in cell_ids:
            self.cluster_labels[self.cluster_labels==cell_id] = self.trash_label
        self.notify_observers()
    
    def merge_cells(self, *cell_ids):
        """merge selected cells. after merging all cells receive the label of the
        first cell"""
        for cell in cell_ids:
            self.cluster_labels[self.cluster_labels==cell]=cell_ids[0]
        self.notify_observers()
        
    def update(self):
        self._cluster(None, self.method, *self.args, **self.kwargs)
        super(ClusterAnalyzer, self).update()
            
    labels = property(read_labels)

class MplPlotComponent(base.Component):
    """Base class for plot components"""
    
    def __init__(self, figsize=(8,6)):
        self.fig = None
        self.figsize = figsize
        super(MplPlotComponent, self).__init__()
        
    def _draw(self, new_figure=False):
        if new_figure or self.fig is None:
            self.fig = plotting.figure(figsize=self.figsize)
            plotting.show()
        self.fig.clf()
        self._plot()
        self.fig.canvas.draw()

    def show(self):
        if not self.fig:
            self._draw()
        #plotting.show()
        
    def update(self):
        if self.fig is not None:
            self._draw()

class ExportCells(base.Component):
    labels_src = base.RequiredFeature("LabelSource", 
                                      base.HasAttributes("labels"))
    marker_src = base.RequiredFeature("SpikeMarkerSource",
                                      base.HasAttributes("events"))
    export_filte = base.RequiredFeature("ExportFilter",
                                        base.HasAttributes("write_spt"))
    
    def export(self):
        labels = labels_src.labels
        spike_idx = marker_src.events
        spt_clust = cluster.cluster2spt(spike_idx, labels)


class PlotFeatures(MplPlotComponent):
    feature_src = base.RequiredFeature("FeatureSource", base.HasAttributes("features"))
    cluster_src = base.RequiredFeature("LabelSource", base.HasAttributes("labels"))
    
    def _plot(self):
        feats = self.feature_src.features
        labels = self.cluster_src.labels
        try:
            plotting.plot_features(feats, labels, fig=self.fig)
        except IndexError:
            pass

       
class PlotSpikes(MplPlotComponent):
    spike_src = base.RequiredFeature("SpikeSource", base.HasAttributes("spikes"))
    cluster_src = base.RequiredFeature("LabelSource", base.HasAttributes("labels"))
    
    def _plot(self):
        spikes = self.spike_src.spikes
        labels = self.cluster_src.labels
        try:
            plotting.plot_spikes(spikes, labels, fig=self.fig)
        except IndexError:
            pass

class Legend(MplPlotComponent):
    cluster_src = base.RequiredFeature("LabelSource", base.HasAttributes("labels"))
    
    def __init__(self):
        super(Legend, self).__init__()
        self.figsize=(1,2)
    
    def _plot(self):
        labels = np.unique(self.cluster_src.labels)
        ax = self.fig.add_axes([0,0,1,1])
        plotting.legend(labels, ax=ax)
    
            