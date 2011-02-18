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
    
    def __init__(self, conf_file, dataset, overwrite=False, f_filter=None):
        BakerlabFilter.__init__(self, conf_file)
        self.dataset = dataset
        self._signal = None
        self._events = None
        self.overwrite = overwrite
        self.f_filter = f_filter
        super(BakerlabSource, self).__init__()
        
    def read_signal(self):
        if self._signal is None:
            self._signal = self.read_sp(self.dataset)
            if self.f_filter is not None:
                filter = sort.extract.Filter("ellip", *self.f_filter)
                self._signal = sort.extract.filter_proxy(self._signal, filter)
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
                 sp_win=(-0.2, 0.8),
                 f_filter=None,
                 align=True):
        self._thresh = thresh
        self.contact = contact
        self.type = type
        self.align = align
        self.resample = resample
        self.sp_win = sp_win
        self.sp_times = None
        self.f_filter = f_filter
        self._est_thresh = None
        super(SpikeDetector, self).__init__()
    
    def _get_threshold(self):
        if self._est_thresh is None:
            return self._thresh
        else:
            return self._est_thresh
        
    def _set_threshold(self, value):
        self._thresh =  value
        self._est_thresh = None
        
    threshold = property(_get_threshold, _set_threshold)
        
    def _detect(self):
        sp = self.waveform_src.signal
        
        if self.f_filter is None:
            filter = None
        else:
            filter = sort.extract.Filter("ellip", *self.f_filter)
            sp = sort.extract.filter_proxy(sp, filter)
        spt = sort.extract.detect_spikes(sp,   edge=self.type,
                                               contact=self.contact,
                                               thresh=self._thresh,
                                               filter=filter)
        self._est_thresh = spt['thresh']
        if self.align:
            self.sp_times = sort.extract.align_spikes(sp, spt, 
                                                      self.sp_win, 
                                                      type=self.type,
                                                      contact=self.contact,  
                                                      resample=self.resample)
        else:
            self.sp_times = spt
    
    def _update(self):
        self._detect()

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
    
    def _update(self):
        self._extract_spikes()
    
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
        self._feature_data = features.combine(feats, norm=self.normalize)
    
    def read_features(self):
        if self._feature_data is None:
            self._calc_features()
        return self._feature_data

    def _update(self):
        self._calc_features()
        super(FeatureExtractor, self)._update()
    
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
        
    def _update(self):
        self._cluster(None, self.method, *self.args, **self.kwargs)
        
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
            self.fig.canvas.mpl_connect("close_event", self._close_callback)
        self.fig.clf()
        self._plot()
        self.fig.canvas.draw()
        
    def _close_callback(self, event):
        self.fig = None

    def show(self):
        if not self.fig:
            self._draw()
        #plotting.show()
        #plotting.show()
        
    def _update(self):
        if self.fig is not None:
            self._draw()

class ExportCells(base.Component):
    labels_src = base.RequiredFeature("LabelSource", 
                                      base.HasAttributes("labels"))
    marker_src = base.RequiredFeature("SpikeMarkerSource",
                                      base.HasAttributes("events"))
    export_filter = base.RequiredFeature("EventsOutput",
                                        base.HasAttributes("events"))
    
    def export(self, mapping=None,overwrite=False):
        labels = self.labels_src.labels
        spike_idx = self.marker_src.events
        self.export_filter.overwrite = overwrite
        export_events = self.export_filter.events
        spt_clust = sort.cluster.split_cells(spike_idx, labels)
        if mapping is None:
            for cell_id, spt in spt_clust.items():
                export_events['cell{0}'.format(cell_id)]=spt
        else:
            for clust_id, cell_id in mapping.items():
                export_events['cell{0}'.format(cell_id)]=spt_clust[clust_id]
                
            


class PlotFeatures(MplPlotComponent):
    feature_src = base.RequiredFeature("FeatureSource", 
                                       base.HasAttributes("features"))
    cluster_src = base.RequiredFeature("LabelSource", 
                                       base.HasAttributes("labels"))
    
    def __init__(self):
        super(PlotFeatures, self).__init__()
        self._showcells = 'all' 
        self._autoscale = False
        
    def _get_autoscale(self):
        return self._autoscale
    
    def _set_autoscale(self, value):
        self._autoscale = value
        if self.fig is not None:
            self._draw()
        
    autoscale = property(_get_autoscale,_set_autoscale, None, 
                         "automatically set plot limits")
    
    def _get_showcells(self):
        return self._showcells
    
    def _set_showcells(self, value):
        if value == 'all':
            self._showcells = value
        else:
            try: 
                it = iter(value)
                self._showcells = value
            except TypeError:
                self._showcells = [value]
        if self.fig is not None:
            self._draw()
        
    show_cells = property(_get_showcells,_set_showcells, None, 
                         "list of labels of cells to plot")
        
    def _get_features(self):
        return self.feature_src.features   
    
    def _update(self):
        self._showcells = 'all'
        super(PlotFeatures, self)._update()

    def _plot(self):
        feats = self._get_features()
        labels = self.cluster_src.labels
        if self.show_cells =='all':
            show_labels = list(np.unique(labels))
            if 0 in show_labels: show_labels.remove(0)
        else:
            show_labels = self.show_cells
        data_range = None if self._autoscale else [0,1]
        plotting.plot_features(feats, labels, show_cells=show_labels, 
                               datarange=data_range, fig=self.fig)
    
class PlotFeaturesTimeline(PlotFeatures):
    spk_time_src =  base.RequiredFeature("SpikeMarkerSource", 
                                         base.HasAttributes("events"))
    
    def _get_features(self):
        spt_dict = self.spk_time_src.events
        feats = self.feature_src.features
        spk_time = sort.features.fetSpTime(spt_dict)
        new_features = sort.features.combine((spk_time,feats), norm=True)
        return new_features
               
class PlotSpikes(MplPlotComponent):
    spike_src = base.RequiredFeature("SpikeSource", base.HasAttributes("spikes"))
    cluster_src = base.RequiredFeature("LabelSource", base.HasAttributes("labels"))
    
    def __init__(self):
        super(PlotSpikes, self).__init__()
        self.show_cells = 'all'
    
    def _update(self):
        self.show_cells = 'all'
        super(PlotSpikes, self)._update()
    
    def _get_showcells(self):
        return self._showcells
    
    def _set_showcells(self, value):
        if value == 'all':
            self._showcells = value
        else:
            try: 
                it = iter(value)
                self._showcells = value
            except TypeError:
                self._showcells = [value]
        if self.fig is not None:
            self._draw()
        
    show_cells = property(_get_showcells,_set_showcells, None, 
                         "list of labels of cells to plot")
        
    def _plot(self):
        spikes = self.spike_src.spikes
        labels = self.cluster_src.labels
        if self.show_cells =='all':
            show_labels = list(np.unique(labels))
            if 0 in show_labels: show_labels.remove(0)
        else:
            show_labels = self.show_cells
        plotting.plot_spikes(spikes, labels, show_cells=show_labels,
                             fig=self.fig)
       
class Legend(MplPlotComponent):
    cluster_src = base.RequiredFeature("LabelSource", base.HasAttributes("labels"))
    
    def __init__(self):
        super(Legend, self).__init__()
        self.figsize=(1,2)
    
    def _plot(self):
        labels = np.unique(self.cluster_src.labels)
        ax = self.fig.add_axes([0,0,1,1])
        plotting.legend(labels, ax=ax)
    
            