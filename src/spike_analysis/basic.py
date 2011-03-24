import matplotlib.pyplot as plt
import numpy as np

def gettrains(spt, stim, win, binsz):
    i = np.searchsorted(stim, spt)
    spt2 = (spt - stim[i - 1])
    bin = np.arange(win[0], win[1], binsz)
    bin = np.concatenate((bin, [bin[-1] + binsz, np.inf]))
    j = np.searchsorted(bin, spt2)
    npoints = len(bin) - 1;
    ntrials = len(stim)
    trains = np.zeros((npoints, ntrials))

    trains[j - 1, i - 1] = 1
    return trains[:-1, :];

def SortSpikes(spt,stim,win=None):
    """Given spike and stimuli times return 2D array with spike trains.
       If win is given only spikes occuring in the time window are 
      returned.
    """
    i = np.searchsorted(stim,spt)
    spt2 = (spt-stim[i-1])
    if win:
        corrected = filter(lambda x: win[1]>x[0]>=win[0], zip(spt2,i))
        spt2 = np.array([x[0] for x in corrected])
        i = np.array([x[1] for x in corrected])
    return [spt2[i==j] for j in xrange(1,len(stim))] 

def plotraster(spt,stim,win=[0,30],ntrials=None,ax=None,height=1.):
    """Creates raster plots of spike trains:
   
         spt - times of spike occurance,
         stim - stimulus times
         win - range of time axis
         ntrials - number of trials to plot (if None plot all)
    """
    if not ntrials: ntrials=len(stim)-1
    if not ax: ax = plt.gca()
 
    spt2 = spt[spt<stim[ntrials]].copy()
    i = np.searchsorted(stim[:ntrials],spt2);
    spt2 = (spt2-stim[i-1]);
  
    plt.vlines(spt2,i,i+height)
    plt.xlim(win)
    plt.ylim((1,ntrials))
    plt.xlabel('time (ms)')
    plt.ylabel('trials')

def plottrains(trains,win=[0,30],ntrials=None,height=1.):
    print "Deprecation: Please use plotRasterTrains insted"
    plotRasterTrains(trains,win,ntrials,height)

def plotRasterTrains(trains,win=[0,30],ntrials=None,height=1.):
    """Creates raster plots of spike trains:
     
         spt - times of spike occurance,
         stim - stimulus times
         win - range of time axis
         ntrials - number of trials to plot (if None plot all)
    """
    if ntrials: 
        trains=trains[:ntrials]
    
    ax = plt.gca();
    ax.set_xlim(win);
    ax.set_ylim((1,ntrials));
    
    lines=[plt.vlines(sp, np.ones(len(sp))*i, np.ones(len(sp))*(i+height)) 
           for i,sp in enumerate(trains) if len(sp)]
    plt.xlabel('time (ms)')
    plt.ylabel('trials')
    
    return lines 

def BinTrains(trains,win,tsamp=0.25):
    """Convert a list of spike trains into a binary sequence"""
    bins = np.arange(win[0], win[1], tsamp)
    trains = [np.histogram(spt,bins)[0][1:] for spt in trains]

    return bins[1:-1], np.array(trains).T

def plotPSTH(spt, stim, win=[0,30], bin=0.25, ax=None, 
             rate=False, **kwargs):
    """Plot peri-stimulus time histogram (PSTH)"""
    i = np.searchsorted(stim+win[0],spt)
    spt2=(spt-stim[i-1])
    bins = np.arange(win[0],win[1],bin)
    (psth,bins) = np.histogram(spt2,bins)
    if not ax:
        ax=plt.gca()
    if rate:
        psth=psth*1./len(stim)/bin*1000.
        ax.set_ylabel('firing rate (Hz)')
    else:
        ax.set_ylabel('number of spikes')
    lines=ax.plot(bins[:-1],psth,**kwargs)
    ax.set_xlabel('time (ms)')
    return lines   

def plotTrainsPSTH(trains,win,bin=0.25, rate=False, **kwargs):
    """Plot peri-stimulus time histogram (PSTH) from a list of spike times
    (trains)"""
    ax = plt.gca()
    time, binned = BinTrains(trains,win,bin)
    psth = np.mean(binned,1)/bin*1000
    plt.plot(time,psth,**kwargs)
    if rate:
        psth=psth*1./len(trains)/bin*1000.
        ax.set_ylabel('firing rate (Hz)')
    else:
        ax.set_ylabel('number of spikes')
    ax.set_xlabel('time (ms)')
    ax.set_ylabel('firing rate (Hz)')

def CalcTrainsPSTH(trains, win, bin=0.25):
    time, binned = BinTrains(trains,win,bin)
    psth = np.mean(binned,1)/bin*1000
    return time, psth

def CalcPSTH(spt,stim,win=[0,30],bin=0.25,ax=None,norm=False,**kwargs):
    """Calculate peri-stimulus time histogram (PSTH).
    Output:
        -- psth - spike counts
        -- bins - bins edges"""
    i = np.searchsorted(stim+win[0],spt)
    spt2 = (spt-stim[i-1])
    bins = np.arange(win[0],win[1],bin)
    (psth, bins) = np.histogram(spt2,bins)
    if norm:
        psth=psth*1000./(len(stim)*bin)
    return psth[1:],bins[1:-1]

def plotPSTHBar(spt,stim,win=[0,30],bin=0.25,**kwargs):
    """Plot peri-stimulus time histogram (PSTH)"""
    i = np.searchsorted(stim+win[0],spt)
    spt2=(spt-stim[i-1])
    bins = np.arange(win[0],win[1],bin)
    (psth, bins) = np.histogram(spt2, bins)
    psth=psth*1./len(stim)/bin*1000.
    plt.gca().bar(bins[:-1], psth, bin, **kwargs)
    plt.xlabel('time (ms)')
    plt.ylabel('firing rate (Hz)')