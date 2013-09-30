from spike_beans import base, components
from nose.tools import ok_, assert_raises
from nose import with_setup
import numpy as np
import json
import os

conf_file = 'test.conf'
el_node = '/Test/s32test01/el1'
cell = 'cell1'

from utils import *

def setup():
    "set up test fixtures"
    base.features = base.FeatureBroker()

def setup_io():
    base.features = base.FeatureBroker()
    file_descr = {"fspike":"{ses_id}{el_id}.sp",
                  "cell":"{ses_id}{el_id}{cell_id}.spt",
                  "dirname":".",
                  "FS":5.E3,
                  "n_contacts":1}
    
    with open(conf_file, 'w') as fp:
         json.dump(file_descr, fp)
         
def teardown_io():
    os.unlink(conf_file)

def teardown():
    "tear down test fixtures"
    
@with_setup(setup, teardown)
def test_renaming():
    class StringSource(base.Component):
        value = "asdf"

    class IntSource(base.Component):
        value = int(5)

    class Consumer(base.Component):
        source = base.RequiredFeature("SourceOne", base.HasAttributes("value"))

    str_src = StringSource()
    int_src = IntSource() 
    consumer = Consumer()

    base.features.Provide("SourceOne", str_src)
    base.features.Provide("SourceTwo", int_src)

    test1 = consumer.source.value == "asdf"
    consumer.source = "SourceTwo" # renaming the requested feature name
    test2 = consumer.source.value == 5

    ok_(test1 and test2)

@with_setup(setup, teardown)
def test_spike_detection():
    base.features.Provide("SignalSource", DummySignalSource())
    detector = components.SpikeDetector(thresh=50.)
    spt = detector.events
    
    source = base.features['SignalSource']
    ok_((np.abs(spt['data']-source.spt)<=1000./source.FS).all())


@with_setup(setup, teardown)
def test_spike_detection_update():
    base.features.Provide("SignalSource", DummySignalSource())
    detector = components.SpikeDetector(thresh=50.)
    spt = detector.events
    detector.threshold = spike_amp*2
    detector.update()
    spt_new = detector.events
    ok_(len(spt_new['data'])==0)

@with_setup(setup_io, teardown_io)
def test_bakerlab_event_read():
    spt_fname = "32test0111.spt"
    spt_data = np.random.randint(0,100, (10,))/200.
    
    (spt_data*200).astype(np.int32).tofile(spt_fname)
    
    src = components.BakerlabSource(conf_file, el_node)
    
    spt_read = src.events[cell]
    
    os.unlink(spt_fname)
    ok_((np.ceil(np.abs(spt_read['data']-spt_data)/200.)<=1).all())
    
@with_setup(setup_io, teardown_io)
def test_bakerlab_event_write():
    spt_data = np.random.randint(0,100, (10,))/200.
    
    conf_file = 'test.conf'
    spt_fname = "32test0111.spt"
    
    src = components.BakerlabSource(conf_file, el_node)
    spt_dict = {"data": spt_data}
    src.events[cell] = spt_dict
    
    ok = os.path.exists(spt_fname)
    os.unlink(spt_fname)
    ok_(ok)
    
@with_setup(setup_io, teardown_io)
def test_export_component():
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    base.features.Provide("LabelSource", DummyLabelSource())
    base.features.Provide("EventsOutput",
                          components.BakerlabSource(conf_file, el_node))
    
    labels = np.unique(base.features['LabelSource'].labels)
    export_comp = components.ExportCells()
    export_comp.export()
    
    for i in labels:
        fname = "32test011{0}.spt".format(i)
        assert os.path.exists(fname)
        os.unlink(fname)
        
@with_setup(setup_io, teardown_io)
def test_export_with_metadata_component():
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    base.features.Provide("LabelSource", DummyLabelSource())
    base.features.Provide("EventsOutput",
                          components.BakerlabSource(conf_file, el_node))
    
    labels = np.unique(base.features['LabelSource'].labels)
    export_comp = components.ExportWithMetadata()
    export_comp.export()
    
    for i in labels:
        fname = "32test011{0}.spt".format(i)
        assert os.path.exists(fname)
        os.unlink(fname)
        
        if i!=0:
            log_fname = "32test011{0}.log".format(i)
            assert os.path.exists(log_fname)
            os.unlink(log_fname)

    
@with_setup(setup_io, teardown_io)
def test_bakerlab_signal_source():
   
    data = np.random.randint(-1000, 1000, (100,))
    fname = "32test011.sp"

    data.astype(np.int16).tofile(fname)
    src = components.BakerlabSource(conf_file, el_node)
    
    sp_read = src.signal
    os.unlink(fname)
    ok_((np.abs(sp_read['data']-data)<=1/200.).all())

@with_setup(setup, teardown)
def test_filter_stack_add_filter_method():
    base.features.Provide("RawSource", DummySignalSource())

    def filter_func(signal):
        new_signal = signal.copy()
        new_signal['data'] += 1.
        return new_signal

    io_filter = components.FilterStack()
    io_filter.add_filter(filter_func)

    ok_((io_filter.signal['data'] > 0.).all())

@with_setup(setup, teardown)
def test_filter_stack_add_filter_string():
    io = DummySignalSource()
    base.features.Provide("RawSource", io)

    io_filter = components.FilterStack()
    io_filter.add_filter("LinearIIR", 800., 300.)
    
    # only tests that filtering has changed a signal
    ok_(np.linalg.norm(io.signal['data'] - io_filter.signal['data']) > 0.)

@with_setup(setup, teardown)
def test_filter_stack_add_filter_attribute_error():
    io = DummySignalSource()
    base.features.Provide("RawSource", io)

    io_filter = components.FilterStack()
    filter_name = "NonExistingFilter"

    assert_raises(AttributeError, io_filter.add_filter, filter_name)

@with_setup(setup, teardown)
def test_filter_stack_add_filter_type_error():
    io = DummySignalSource()
    base.features.Provide("RawSource", io)

    io_filter = components.FilterStack()
    filter_argument = None

    assert_raises(TypeError, io_filter.add_filter, filter_argument)

@with_setup(setup, teardown)
def test_spike_extractor():
    base.features.Provide("SignalSource", DummySignalSource())
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    
    sp_waves = components.SpikeExtractor().spikes
    mean_wave = sp_waves['data'][:,:,0].mean(1) 
    time = sp_waves['time']
    true_spike = spike_amp*((time>=0) & (time<spike_dur))
    ok_(np.sum(np.abs(mean_wave-true_spike))<0.01*spike_amp)

@with_setup(setup, teardown)
def test_feature_extractor():
    base.features.Provide("SignalSource",      DummySignalSource())
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    base.features.Provide("SpikeSource",       components.SpikeExtractor())
    
    feat_comp = components.FeatureExtractor(normalize=False)
    feat_comp.add_feature("P2P")
    features = feat_comp.features
    
    ok_((features['data']==spike_amp).all())       

@with_setup(setup, teardown)
def test_feature_extractor_hide_features():
    base.features.Provide("SignalSource",      DummySignalSource())
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    base.features.Provide("SpikeSource",       components.SpikeExtractor())
    
    feat_comp = components.FeatureExtractor(normalize=False)
    feat_comp.add_feature("P2P")
    feat_comp.add_feature("SpIdx")
    feat_comp.hide_features("Sp*")
    feat_comp.update()
    features = feat_comp.features

    test1 = features['names'] == ["P2P:Ch0:P2P"]  # it's P2P
    test2 = features['data'].shape[1] == 1
    test3 = (features['data'] == spike_amp).all()  # with corresponding data
    
    ok_(test1 and test2 and test3)       

@with_setup(setup, teardown)
def test_feature_extractor_hide_features_not_found():
    base.features.Provide("SignalSource",      DummySignalSource())
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    base.features.Provide("SpikeSource",       components.SpikeExtractor())
    
    feat_comp = components.FeatureExtractor(normalize=False)
    feat_comp.add_feature("P2P")
    feat_comp.update()

    assert_raises(ValueError, feat_comp.hide_features, "NonExistingFeature")

@with_setup(setup, teardown)
def test_feature_extractor_unhide_features():
    base.features.Provide("SignalSource",      DummySignalSource())
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    base.features.Provide("SpikeSource",       components.SpikeExtractor())
    
    feat_comp = components.FeatureExtractor(normalize=False)
    feat_comp.add_feature("P2P")
    feat_comp.add_feature("SpIdx")
    feat_comp.hide_features("*")
    feat_comp.update()
    feat_comp.unhide_features("*Ch0*")
    feat_comp.update()
    features = feat_comp.features

    test1 = features['names'] == ["P2P:Ch0:P2P"]   # it's PC1
    test2 = features['data'].shape[1] == 1  # with corresponding data
    
    ok_(test1 and test2)       

@with_setup(setup, teardown)
def test_feature_extractor_unhide_features_not_found():
    base.features.Provide("SignalSource",      DummySignalSource())
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    base.features.Provide("SpikeSource",       components.SpikeExtractor())
    
    feat_comp = components.FeatureExtractor(normalize=False)
    feat_comp.add_feature("P2P")
    feat_comp.add_feature("SpIdx")
    feat_comp.hide_features("*SpIdx")
    feat_comp.update()

    assert_raises(ValueError, 
            feat_comp.unhide_features,
            "NonExistingFeature")

@with_setup(setup, teardown)
def test_cluster_component():
    base.features.Provide("FeatureSource", DummyFeatureExtractor())
    
    cluster_comp = components.ClusterAnalyzer("k_means", 2)
    labels = cluster_comp.read_labels()
    
    ok = (((labels[:n_spikes]==1).all() & (labels[n_spikes:]==2).all()) |
         ((labels[:n_spikes]==2).all() & (labels[n_spikes:]==1).all()))
    ok_(ok)
    

@with_setup(setup, teardown)
def test_cluster_component_relabel():
    base.features.Provide("FeatureSource", RandomFeatures())
    
    cluster_comp = components.ClusterAnalyzer("k_means", 5)
    labs = cluster_comp.labels
    cluster_comp.delete_cells(1,2,3,4)
    cluster_comp.relabel()
    
    labels = np.unique(cluster_comp.labels)
    labels.sort()
    
    ok_((labels==np.array([0,1])).all())

@with_setup(setup, teardown)
def test_cluster_component_smart_update():
    feature_comp = base.register("FeatureSource", DummyFeatureExtractor())
    cluster_comp = components.ClusterAnalyzer("k_means", 2)
    labs = cluster_comp.labels # this is a workaround to call _cluster()

    # cluster_comp should NOT recluster when updated, if the number of
    # spikes didn't change
    cluster_comp.delete_cells(1) # modify lablels
    labels_orig = cluster_comp.labels.copy()

    feature_comp.add_feature('new_feature')
    feature_comp.update()
    labels_new_feat = cluster_comp.labels.copy()

    # cluster_comp should recluster when updated, if the number of
    # spikes DID change
    feature_comp.add_spikes(10)
    feature_comp.update()
    labels_new_spikes = cluster_comp.labels.copy()

    test1 = (labels_orig == labels_new_feat).all()
    test2 = len(labels_new_spikes) != len(labels_new_feat)

    ok_(test1 and test2)

@with_setup(setup, teardown)
def test_cluster_component_methods_before_labels_requested():
    # issue #75
    # The followig methods should not fail when called before lebels are
    # externally requested (cluster_labels = None)

    base.features.Provide("FeatureSource", RandomFeatures())
    cluster_comp = components.ClusterAnalyzer("k_means", 2)

    cluster_comp.delete_cells(1)

    cluster_comp.cluster_labels = None
    cluster_comp.delete_spikes([0])

    cluster_comp.cluster_labels = None
    cluster_comp.merge_cells(1, 2)

    cluster_comp.cluster_labels = None
    cluster_comp.relabel()

    ok_(cluster_comp.cluster_labels is not None)

@with_setup(setup, teardown)
def test_truncated_spikes_from_end():
    signal_src = DummySignalSource()
    signal_src._spikes = signal_src._spikes[:, :-period/1000.*FS*2.5]
    base.features.Provide("SignalSource",      signal_src)
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    sp_waves = components.SpikeExtractor().spikes

    correct_mask = np.ones(n_spikes-2).astype(np.bool)
    correct_mask[-1] = False
    ok_((sp_waves['is_valid'] == correct_mask).all())    
   
@with_setup(setup, teardown)
def test_truncated_spikes_from_begin():
    signal_src = DummySignalSource()
    detector = DummySpikeDetector()
    spt = detector._spt_data['data']
    detector._spt_data['data'] = np.insert(spt, 0, 0)
    base.features.Provide("SignalSource",      signal_src)
    base.features.Provide("SpikeMarkerSource", detector)
    sp_waves = components.SpikeExtractor().spikes
    
    correct_mask = np.ones(n_spikes-1).astype(np.bool)
    correct_mask[0] = False
    ok_((sp_waves['is_valid'] == correct_mask).all())

@with_setup(setup, teardown)    
def test_null_labels_returned_for_truncated_spikes():
    signal_src = DummySignalSource()
    signal_src._spikes = signal_src._spikes[:, :-period/1000.*FS*2.5]
    
    base.features.Provide("SignalSource",      signal_src)
    base.features.Provide("SpikeMarkerSource", DummySpikeDetector())
    base.features.Provide("SpikeSource",       components.SpikeExtractor())
    base.features.Provide("FeatureSource",     components.FeatureExtractor(normalize=False))
    base.features.Provide("ClusterAnalyzer",   components.ClusterAnalyzer("k_means", 1))
    base.features['FeatureSource'].add_feature("P2P")
    cl = base.features["ClusterAnalyzer"].labels
    
    true_labels = np.ones(n_spikes-2)
    true_labels[-1]= 0
    
    ok_((cl==true_labels).all())

@with_setup(setup, teardown)
def test_propagate_truncate_to_features():
    
    spike_src = DummySpikeSource()
    sp_waves = spike_src._sp_waves
    is_valid = np.ones(sp_waves['data'].shape[1]).astype(bool)
    spike_src._sp_waves['is_valid'] = is_valid

    base.features.Provide("SpikeSource",  spike_src)
    
    feat_comp = components.FeatureExtractor(normalize=False)
    feat_comp.add_feature("P2P")
    
    features = feat_comp.features
    
    ok_((features['is_valid']==is_valid).all())  
    

@with_setup(setup, teardown)
def test_pipeline_update():
    base.features.Provide("SignalSource",      DummySignalSource())
    base.features.Provide("SpikeMarkerSource", components.SpikeDetector(thresh=spike_amp/2.))
    base.features.Provide("SpikeSource",       components.SpikeExtractor())
    base.features.Provide("FeatureSource",     components.FeatureExtractor(normalize=False))
    base.features.Provide("ClusterAnalyzer",   components.ClusterAnalyzer("k_means", 2))
    
    base.features['FeatureSource'].add_feature("P2P")
    
    cl1 = base.features["ClusterAnalyzer"].labels
    base.features["SignalSource"].update()
    cl2 = base.features["ClusterAnalyzer"].labels
    ok_(len(cl1)==98)
    ok_(len(cl2)==48)
    
