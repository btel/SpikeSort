'''
Created on Feb 9, 2011

@author: bartosz
'''

import spike_sort as sort

import base

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