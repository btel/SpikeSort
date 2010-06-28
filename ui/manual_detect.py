#!/usr/bin/env python
#coding=utf-8

import numpy as np
import matplotlib.pyplot as plt
from ..core import extract
import patterns

i = 0

def find_spikes(spike_dict, spt_dict, trigger_dict,  win, sp_win=[-0.2, 0.8]):

    def next_plot(inc=1):
        global i
        ax.cla()
        i=(i+inc) % sp_waves.shape[1]
        ax.plot(sp_time, sp_waves[:,i], picker=5)
        if len(trains_new[i])>0:
            ax.plot(trains_new[i],np.ones(len(trains_new[i]))*0.8,
                    "*")
        if i in missed_trains:
            ax.plot(missed_trains[i],np.ones(len(missed_trains[i]))*0.8,
                    "r*")
        ax.set_ylim(0,1)
        ax.set_title("%d/%d" % (i+1, sp_waves.shape[1]))
        fig.canvas.draw()

    def onclick(event):
        if event.button == 1:
            if event.inaxes is None:
                next_plot()
        elif event.button == 3:
            next_plot(-1)


    def onpick(event):
        if event.mouseevent.button == 1:
            #append a new spike
            _win = (np.array(sp_win)/1000.*FS).astype(int)
            sp_idx = (stim[i] + win[0])/1000.*FS + event.ind[0]
            spike = sp_raw[int(sp_idx)+_win[0]:int(sp_idx)+_win[1]]
        
            missed_spikes.append(sp_idx*1000./FS)
            if not i in missed_trains:
                missed_trains[i] = []
            missed_trains[i].append(sp_idx*1000./FS-stim[i])
            ax2.plot(spike, 'k')
            fig2.canvas.draw()
            next_plot()
        elif event.mouseevent.button == 3:
            #implement deleting a spike
            #for now just go to the previous trial
            next_plot(-1)



    stim = trigger_dict['data']
    spt = spt_dict['data']
    sp_raw = spike_dict['data']
    FS = spike_dict['FS']
   
    trains_new = patterns.SortSpikes(spt, stim, win)
   
    spike_trials = extract.extract_spikes(spike_dict, trigger_dict, win)
    sp_waves = spike_trials['data']
    sp_time = spike_trials['time']
    sp_waves = (sp_waves - sp_waves.min())/(sp_waves.max()-sp_waves.min())

    missed_spikes = []
    missed_trains = {}
    
    fig = plt.figure()
    ax = plt.subplot(111)
    next_plot()
    fig.canvas.mpl_connect("button_press_event", onclick)
    fig.canvas.mpl_connect('pick_event', onpick)

    fig2 = plt.figure()
    ax2 = plt.subplot(111)
    plt.show()

    missed_spikes = np.array(missed_spikes)

    return {"data":missed_spikes}


