'''
Created on Feb 9, 2011

@author: bartosz
'''

import spike_sort as sort

import base
from spike_sort.io.filters import BakerlabFilter, PyTablesFilter
from spike_sort import features, filters
from spike_sort.ui import plotting
from spike_sort.ui import zoomer
from spike_analysis import dashboard
import numpy as np

from collections import OrderedDict
import fnmatch
import re


class GenericSource(base.Component):
    def __init__(self, dataset, overwrite=False):
        self.dataset = dataset
        self._signal = None
        self._events = None
        self.overwrite = overwrite
        super(GenericSource, self).__init__()

    def read_signal(self):
        if self._signal is None:
            self._signal = self.read_sp(self.dataset)
        return self._signal

    def read_events(self, cell):
        node = '/'.join((self.dataset, cell))
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


class BakerlabSource(GenericSource, BakerlabFilter):
    def __init__(self, conf_file, dataset, overwrite=False):
        GenericSource.__init__(self, dataset, overwrite)
        BakerlabFilter.__init__(self, conf_file)


class PyTablesSource(GenericSource, PyTablesFilter):
    def __init__(self, h5file, dataset, overwrite=False):
        GenericSource.__init__(self, dataset, overwrite)
        PyTablesFilter.__init__(self, h5file)


class FilterStack(base.Component):
    raw_src = base.RequiredFeature("RawSource",
                                   base.HasAttributes("signal"))

    def __init__(self):
        self._filters = []
        self._signal = None
        super(FilterStack, self).__init__()

    def add_filter(self, filt, *args, **kwargs):
        # type checking
        if hasattr(filt, "__call__"):
            self._filters.append(lambda signal: filt(signal, *args, **kwargs))
        elif isinstance(filt, str):
            try:
                filter_func = filters.__getattribute__("flt" + filt)
                self._filters.append(
                    lambda signal: filter_func(signal, *args, **kwargs))
            except AttributeError:
                raise AttributeError("No such method found in 'core.filters':"
                                     + " flt" + filt)
        else:
            raise TypeError(("Unsupported argument type: %s." % type(filt)) +
                            "Only string or callable are accepted")

    def read_signal(self):
        if self._signal is None:
            self._signal = self.raw_src.signal
            for filt in self._filters:
                self._signal = filt(self._signal)
        return self._signal

    signal = property(read_signal)


class SpikeDetector(base.Component):
    """Detect Spikes with alignment"""
    waveform_src = base.RequiredFeature("SignalSource",
                                        base.HasAttributes("signal"))

    def __init__(self, thresh='auto',
                 contact=0,
                 type='max',
                 resample=1,
                 sp_win=(-0.2, 0.8),
                 align=True):
        self._thresh = thresh
        self.contact = contact
        self.type = type
        self.align = align
        self.resample = resample
        self.sp_win = sp_win
        self.sp_times = None
        self._est_thresh = None
        super(SpikeDetector, self).__init__()

    def _get_threshold(self):
        return self._est_thresh or self._thresh

    def _set_threshold(self, value):
        self._thresh = value
        self._est_thresh = None

    threshold = property(_get_threshold, _set_threshold)

    def _detect(self):
        sp = self.waveform_src.signal

        spt = sort.extract.detect_spikes(sp, edge=self.type,
                                         contact=self.contact,
                                         thresh=self._thresh)
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

    def __init__(self, sp_win=[-0.2, 0.8]):
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
        self.feature_methods = OrderedDict()
        self._feature_data = None
        self._hidden_features = []
        self.normalize = normalize
        super(FeatureExtractor, self).__init__()

    def add_feature(self, name, *args, **kwargs):
        func_name = "fet" + name
        _func = features.__getattribute__(func_name)
        func = lambda x: _func(x, *args, **kwargs)
        name = features.create_method_name(name, self.feature_methods.keys())

        self.feature_methods[name] = func

    def hide_features(self, pattern):
        """Hide featues, with names matching `pattern`.

        Patterns are Unix shell style:

        *       matches everything
        ?       matches any single character
        [seq]   matches any character in seq
        [!seq]  matches any char not in seq

        Parameters
        ----------
        pattern : string
            search pattern
        """
        if not isinstance(pattern, basestring):
            raise TypeError("Wrong pattern type: %s. Must be string"
                    % type(pattern))

        fts_to_delete = fnmatch.filter(self.features['names'], pattern)

        if not fts_to_delete:
            raise ValueError("No matching features found")
        
        self._hidden_features = list(
                set(fts_to_delete) | set(self._hidden_features))

    def unhide_features(self, pattern):
        """Recover features matching `pattern`, if they were previously
        hidden with `hide_features`

        Parameters
        ----------
        pattern : string
            search pattern
        """
        fts_to_undelete = fnmatch.filter(self._hidden_features, pattern)
        if not fts_to_undelete:
            raise ValueError(
                "Error: features matching '%s' were never hidden. \
                 Use `add_features` to add new features" % pattern)

        for ft in fts_to_undelete:
            self._hidden_features.remove(ft)

    def clear_selection(self):
        """ Bring back all previously hidden features
        """
        self._hidden_features = []

    def _calc_features(self):
        spikes = self.spikes_src.spikes
        feats = [f(spikes) for f in self.feature_methods.values()]
        ft_data = features.combine(feats,
                norm=self.normalize,
                feat_method_names = self.feature_methods.keys())
        
        # Filter feature_data to remove _hidden_features.
        # This routine is O(n^2), to deal woth possible repetitions
        names, idx = [], []
        for i, name in enumerate(ft_data['names']):
            if name not in self._hidden_features:
                names.append(name)
                idx.append(i)

        ft_data['names'] = names
        ft_data['data'] = ft_data['data'][:, idx]

        self._feature_data = ft_data

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
        self.force_recluster_on_update = False
        super(ClusterAnalyzer, self).__init__()
        self.use_features = 'all'

    def _cluster(self, idx, method, *args, **kwargs):
        feature_data = self.feature_src.features
        use_features = self.use_features

        if use_features != 'all' and use_features is not None:
            names = list(feature_data['names'])
            try:
                ii = np.array([names.index(l) for l in use_features])
                feature_data = feature_data.copy()
                feature_data['data'] = feature_data['data'][:, ii]
                feature_data['names'] = feature_data['names'][ii]
            except ValueError:
                raise ValueError("Feature {0} does not exist" % l)

        if idx is not None:
            new_features = sort.features.select_spikes(feature_data, idx)
            clust_idx = sort.cluster.cluster(method, new_features, *args,
                                             **kwargs)
            all_labels = set(range(1, 100))
            used_labels = set(np.unique(self.cluster_labels))
            free_labels = list(all_labels - used_labels)
            free_labels.sort(reverse=True)
            new_clust_idx = np.zeros(len(clust_idx), dtype="int16")
            for l in np.unique(clust_idx):
                new_label = free_labels.pop()
                new_clust_idx[clust_idx == l] = new_label

            self.cluster_labels[idx] = new_clust_idx
        else:
            clust_idx = sort.cluster.cluster(method, feature_data, *args,
                                             **kwargs)
            self.cluster_labels = clust_idx + 1

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
            self.cluster_labels[self.cluster_labels == l] = i + 1
        self.notify_observers()

    def recluster(self, label=None, method=None, *args, **kwargs):
        """repeat clustering of a selected cell"""
        if method is None:
            method = self.method
        if not args:
            args = self.args
        if not kwargs:
            kwargs = self.kwargs
        if label is None:
            idx = None
        else:
            idx = self.labels == label
        self._cluster(idx, method, *args, **kwargs)
        self.notify_observers()

    def delete_cells(self, *cell_ids):
        """move selected labels to thrash (cluster 0).
        if 'all' thrash all cells """
        if len(cell_ids) == 1 and cell_ids[0] == 'all':
            cell_id = np.unique(self.labels)
        for cell_id in cell_ids:
            idx = self.labels == cell_id
            self.cluster_labels[idx] = self.trash_label
        self.notify_observers()

    def delete_spikes(self, idx_list):
        """moves specified spikes to trash

        Parameters
        ----------
        idx_list : list
            list of spike indices to remove

        """
        self.read_labels()  # prevent NoneType assignment
        self.cluster_labels[idx_list] = self.trash_label
        self.notify_observers()

    def merge_cells(self, *cell_ids):
        """merge selected cells. after merging all cells receive the
        label of the first cell"""
        for cell in cell_ids:
            idx = self.labels == cell
            self.cluster_labels[idx] = cell_ids[0]
        self.notify_observers()

    def _update(self):
        n_spikes_old = self.labels.shape[0]
        n_spikes_new = self.feature_src.features['data'].shape[0]

        if (n_spikes_new != n_spikes_old) or self.force_recluster_on_update:
            self._cluster(None, self.method, *self.args, **self.kwargs)

    labels = property(read_labels)


class MplPlotComponent(base.Component):
    """Base class for plot components"""

    def __init__(self, figsize=(8, 6)):
        self.fig = None
        self.figsize = figsize
        super(MplPlotComponent, self).__init__()

    def _draw(self, new_figure=False):
        if new_figure or self.fig is None:
            self.fig = plotting.figure(figsize=self.figsize)
            self.fig.canvas.mpl_connect("close_event", self._close_callback)
            self.zoomer = zoomer.Zoomer(plotting.plt, self.fig)
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


class PlotFeatures(MplPlotComponent):
    feature_src = base.RequiredFeature("FeatureSource",
                                       base.HasAttributes("features"))
    cluster_src = base.OptionalFeature("LabelSource",
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

    autoscale = property(_get_autoscale, _set_autoscale, None,
                         "automatically set plot limits")

    def _get_showcells(self):
        return self._showcells

    def _set_showcells(self, value):
        if value == 'all':
            self._showcells = value
        else:
            try:
                # it = iter(value)
                self._showcells = value
            except TypeError:
                self._showcells = [value]
        if self.fig is not None:
            self._draw()

    show_cells = property(_get_showcells, _set_showcells, None,
                          "list of labels of cells to plot")

    def _get_features(self):
        return self.feature_src.features

    def _update(self):
        self._showcells = 'all'
        super(PlotFeatures, self)._update()

    def _plot(self):
        feats = self._get_features()

        if self.cluster_src is not None:
            labels = self.cluster_src.labels
            if self.show_cells == 'all':
                show_labels = list(np.unique(labels))
                if 0 in show_labels:
                    show_labels.remove(0)
            else:
                show_labels = self.show_cells
        else:
            labels = None
            show_labels = None

        data_range = None if self._autoscale else [0, 1]
        plotting.plot_features(feats, labels, show_cells=show_labels,
                               datarange=data_range, fig=self.fig)

from spike_sort.ui import spike_browser


class SpikeBrowser(base.Component):
    raw_src = base.RequiredFeature("SignalSource")
    spt_src = base.RequiredFeature("SpikeMarkerSource")
    label_src = base.OptionalFeature("LabelSource",
                                     base.HasAttributes("labels"))

    def __init__(self):
        super(SpikeBrowser, self).__init__()
        self.win = 50
        self.frame = None
        self._showcells = 'all'

    def _on_close(self):
        self.frame.root.destroy()
        self.frame = None

    def _set_data(self):
        if self.label_src is None:
            self._set_data_without_labels()
        else:
            self._set_data_with_labels()

    def _set_data_without_labels(self):
        sp_data = self.raw_src.signal
        spike_time = self.spt_src.events
        self.browser.winsz = self.win
        self.browser.set_data(sp_data)

        self.browser.set_spiketimes(spike_time)

    def _set_data_with_labels(self):
        sp_data = self.raw_src.signal
        spike_time = self.spt_src.events
        labels = self.label_src.labels
        self.browser.winsz = self.win
        self.browser.set_data(sp_data)
        if self._showcells == 'all':
            self.browser.set_spiketimes(spike_time, labels)
        else:
            i = np.in1d(labels, self._showcells)
            spike_time = spike_time.copy()
            spike_time['data'] = spike_time['data'][i]
            self.browser.set_spiketimes(spike_time, labels[i],
                                        np.unique(labels))

    def _get_showcells(self):
        return self._showcells

    def _set_showcells(self, value):
        if value == 'all':
            self._showcells = value
        else:
            try:
                self._showcells = value
            except TypeError:
                self._showcells = [value]
        self.update()

    show_cells = property(_get_showcells, _set_showcells, None,
                          "list of labels of cells to plot")

    def _draw(self):
        self.frame = spike_browser.PlotWithScrollBarTk()
        self.frame.root.protocol("WM_DELETE_WINDOW", self._on_close)
        self.browser = spike_browser.SpikeBrowserUI(self.frame)
        self._set_data()

    def _update(self):
        if self.frame:
            self._set_data()
            self.browser.draw_plot()

    def show(self):
        if not self.frame:
            self._draw()


class PlotFeaturesTimeline(PlotFeatures):
    spk_time_src = base.RequiredFeature("SpikeMarkerSource",
                                        base.HasAttributes("events"))

    def _get_features(self):
        spt_dict = self.spk_time_src.events
        feats = self.feature_src.features
        spk_time = sort.features.fetSpTime(spt_dict)
        new_features = sort.features.combine((spk_time, feats), norm=True)
        return new_features


class PlotSpikes(MplPlotComponent):
    spike_src = base.RequiredFeature("SpikeSource",
                                     base.HasAttributes("spikes"))
    cluster_src = base.OptionalFeature("LabelSource",
                                       base.HasAttributes("labels"))

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
                self._showcells = value
            except TypeError:
                self._showcells = [value]
        if self.fig is not None:
            self._draw()

    show_cells = property(_get_showcells, _set_showcells, None,
                          "list of labels of cells to plot")

    def _plot(self):
        spikes = self.spike_src.spikes

        if self.cluster_src:
            labels = self.cluster_src.labels
            if self.show_cells == 'all':
                show_labels = list(np.unique(labels))
                if 0 in show_labels:
                    show_labels.remove(0)
            else:
                show_labels = self.show_cells
        else:
            labels = None
            show_labels = None

        plotting.plot_spikes(spikes, labels, show_cells=show_labels,
                             fig=self.fig)


class Legend(MplPlotComponent):
    cluster_src = base.RequiredFeature("LabelSource",
                                       base.HasAttributes("labels"))

    def __init__(self):
        super(Legend, self).__init__()
        self.figsize = (2, 2)

    def _plot(self):
        labels = np.unique(self.cluster_src.labels)
        ax = self.fig.add_axes([0, 0, 1, 1])
        plotting.legend(labels, ax=ax)

    def _draw(self, *args, **kwargs):
        # hide the toolbar
        tb = plotting.plt.rcParams['toolbar']
        plotting.plt.rcParams['toolbar'] = 'None'
        super(Legend, self)._draw(*args, **kwargs)
        plotting.plt.rcParams['toolbar'] = tb


class ExportCells(base.Component):
    labels_src = base.RequiredFeature("LabelSource",
                                      base.HasAttributes("labels"))
    marker_src = base.RequiredFeature("SpikeMarkerSource",
                                      base.HasAttributes("events"))
    export_filter = base.RequiredFeature("EventsOutput",
                                         base.HasAttributes("events"))

    def export(self, mapping=None, overwrite=False, metadata='default'):
        labels = self.labels_src.labels
        spike_idx = self.marker_src.events
        self.export_filter.overwrite = overwrite
        export_events = self.export_filter.events
        spt_clust = sort.cluster.split_cells(spike_idx, labels)

        if metadata == 'default':
            md = self.get_metadata()
        else:
            md = metadata

        if mapping is None:
            for cell_id, spt in spt_clust.items():
                if md and cell_id != 0:
                    spt['metadata'] = md
                export_events['cell{0}'.format(cell_id)] = spt
        else:
            for clust_id, cell_id in mapping.items():
                spt = spt_clust[clust_id]
                if md and cell_id != 0:
                    spt['metadata'] = md
                export_events['cell{0}'.format(cell_id)] = spt

    def get_metadata(self):
        return {}


class ExportWithMetadata(ExportCells):
    marker_src = base.RequiredFeature("SpikeMarkerSource",
                                      base.HasAttributes("events",
                                                         "threshold",
                                                         "type",
                                                         "sp_win"))

    def get_metadata(self):
        metadata = {'contact': self.marker_src.contact,
                    'threshold': self.marker_src.threshold,
                    'type': self.marker_src.type,
                    'sp_win': self.marker_src.sp_win}
        return metadata


class Dashboard(MplPlotComponent):
    labels_src = base.RequiredFeature("LabelSource",
                                      base.HasAttributes("labels"))
    marker_src = base.RequiredFeature("SpikeMarkerSource",
                                      base.HasAttributes("events"))
    export_filter = base.RequiredFeature("EventsOutput",
                                         base.HasAttributes("dataset"))

    def _plot(self):
        stim = self.export_filter.read_spt("/".join(
                                           (self.export_filter.dataset,
                                            'stim')))

        labels = self.labels_src.labels
        spike_idx = self.marker_src.events

        spt_all = sort.cluster.split_cells(spike_idx, labels)

        if self.cell not in spt_all.keys():
            old_cell = self.cell
            for c in range(max(self.labels_src.labels) + 1)[::-1]:
                if c in spt_all:
                    self.cell = c
                    print("Dashboard: cell %s doesn't exist, plotting cell %s"
                          % (old_cell, self.cell))
                    break

        dataset = {'spt': spt_all[self.cell]['data'],
                   'stim': stim['data'],
                   'ev': []}
        dashboard.plot_dataset(dataset, self.fig)

    def show(self, cell):
        self.cell = cell
        self._draw()
