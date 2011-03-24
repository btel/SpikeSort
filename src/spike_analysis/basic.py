from pylab import *
from scipy import *
import pdb
import unittest
import numpy as np

def BurstIntervals(isi,maxISI):
    z=hstack((False,isi<=maxISI,False))
    (i,)=where(~z)
    blens=hstack(([0],diff(i)))
    z2=cumsum(~z)
    z2[i]=0
    seq=blens[z2]
    return seq[1:-1]

def FilterSpikes(spt,stim,win):
     i=searchsorted(stim,spt);
     spt2=(spt-stim[i-1]);
     spikes=spt[(spt2<win[1]) & (spt2>=win[0])]
     return spikes

def ShuffleSpikes(spt,stim):
    """Perform exchange resampling"""
    i=searchsorted(stim,spt);
    sh_spt=spt.copy()
    sh_spt=(sh_spt-stim[i-1]);
    random.shuffle(i)
    sh_spt=sh_spt+stim[i-1]
    sh_spt.sort()
    return sh_spt

def EventBorders(trains,bin=0.2,kernel_width=0.25,N=30,threshold=0.15):
    """Segementation of PSTH histograms"""
    spt=concatenate(trains)
    win=(min(spt),max(spt))
    bins=arange(win[0],win[1]+bin,bin)
    (psth,bins)=histogram(spt,bins-bin/2)
    
    t=arange(N)*bin
    y=exp(-(t-N/2*bin)**2/(2*kernel_width**2));
    psthS=convolve(psth,y/sum(y))[N/2:];

    D_psthS=diff(psthS)
    extrema,=where(D_psthS[:-1]*D_psthS[1:]<0)
    exval=psthS[extrema+1]
    base=exval[1:-1]
    ratio=sqrt((exval[:-2]-base)*(exval[2:]-base))/sqrt(exval[:-2]*exval[2:])
    k,=where((ratio>=threshold) & (base<exval[:-2]))
    events=(extrema[k+1]+1)*bin+win[0]
    events=r_[win[0],events,win[1]]
   
    return events,psth

def FindClasses(trains,ev):
    """Classify spike trains"""
    tr_cl =[x[(x>ev[0]) & (x<ev[-1]) ] for x in trains]  
    cl=array([sum(2**(unique(searchsorted(ev,y))-1)) for y in tr_cl])
    cl[cl<0]=0
    return cl

def CountClasses(cl, lab=None):
    """count the number of class instances"""

    if lab is None:
        lab = np.unique(cl)
    n = np.array([np.sum(cl==l) for l in lab])

    return n, lab


def EventInformation(cl,event):
    bins=arange(cl.min(),cl.max()+1)
    p,binS=histogram(cl,bins)
    p=p*1./sum(p)
    b= (bins & event)==event
    p_x=(p[b]+p[~b])
    p_y=array([sum(p[b]),sum(p[~b])])
    
    h_xy=-(p*log2(p))
    h_x=-(p_x*log2(p_x))
    h_y=-(p_y*log2(p_y))
    
    m_xy=(sum(h_x[~isnan(h_x)])+ 
         sum(h_y[~isnan(h_y)])- 
         sum(h_xy[~isnan(h_xy)]))
    return m_xy

def gettrains(spt,stim,win,binsz):
  i=searchsorted(stim,spt)
  spt2=(spt-stim[i-1])
  #spt2=compress(logical_and(spt2>win[0],spt2<win[1]),spt2)
  bin=arange(win[0],win[1],binsz)
  bin=concatenate((bin,[bin[-1]+binsz,inf]))
  j=searchsorted(bin,spt2)
  npoints=len(bin)-1;
  ntrials=len(stim)
  trains=zeros((npoints,ntrials));
  #for i in stim:
  #  ind=(compress(logical_and(spt>stim[i]+win[0]*200, spt<stim[i]+win[1]*200), spt)-stim[i])/10;
  #  trains[i,ind]=1;
  trains[j-1,i-1]=1
  return trains[:-1,:];

def SortSpikes(spt,stim,win=None):
  """Given spike and stimuli times return 2D array with spike trains.
     If win is given only spikes occuring in the time window are 
     returned.
  """
  i=searchsorted(stim,spt);
  spt2=(spt-stim[i-1]);
  if win:
     corrected=filter(lambda x: win[1]>x[0]>=win[0], zip(spt2,i))
     spt2=array([x[0] for x in corrected])
     i=array([x[1] for x in corrected])
  return [spt2[i==j] for j in xrange(1,len(stim))] 

def plotraster(spt,stim,win=[0,30],ntrials=None,ax=None,height=1.):
  """Creates raster plots of spike trains:
   
       spt - times of spike occurance,
       stim - stimulus times
       win - range of time axis
       ntrials - number of trials to plot (if None plot all)
  """
  if not ntrials: ntrials=len(stim)-1
  if not ax: ax = gca()
 
  spt2=spt[spt<stim[ntrials]].copy()
  i=searchsorted(stim[:ntrials],spt2);
  spt2=(spt2-stim[i-1]);
  
  vlines(spt2,i,i+height)
  xlim(win)
  ylim((1,ntrials))
  xlabel('time (ms)')
  ylabel('trials')

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
  
  ax=gca();
  ax.set_xlim(win);
  ax.set_ylim((1,ntrials));

  lines=[vlines(sp,ones(len(sp))*i,ones(len(sp))*(i+height)) for i,sp in enumerate(trains) if len(sp)]
  xlabel('time (ms)')
  ylabel('trials')

  return lines 

def gausswin(n,sigma):
   x=arange(-n/2,n/2)/20.0;
   y=1/(sqrt(2*pi)*sigma)*exp(-x**2/sigma**2);
   return y;

def sigmoid(x,m,tau):
  return 1/(1+exp(-(x-m)/tau));

def residuals(p, y, x):  
   m,tau = p  
   err = y-sigmoid(x,m,tau)  
   return err  

def ReorderMatrix(x,classes,dist=[]):
  """Reorder similarity matrix so that trains belonging to one class are 
  near each other.
  If distances from the centres of the cluster (dist) are also given, it
  additionaly sorts the trials in the order of increasing distance. """
 
  if len(dist)>0:
    i=lexsort([dist,classes])    
  else:
    i = classes.argsort()
  x2=x.copy()
  x2=x2[:,i]
  x2=x2[i,:]
  return x2

def ReorderClassesTrains(trains,classes,dist=[]):
  """Reorder raw (unsorted) spike trains according to which class they
     belong. 
     Use the function to plot sorted raster plot (plotraster)
  """
  
  if len(dist)>0:
    k=lexsort([dist,classes])    
  else:
    k = classes.argsort();
  
  trains = [trains[i] for i in k]
  return trains

def ReorderTrains(spt,stim,classes,dist=[]):
    raise "Deprecation: Please use ReordedClasses"
    

def ReorderClasses(spt,stim,classes,dist=[]):
  """Reorder raw (unsorted) spike trains according to which class they
     belong. 
     Use the function to plot sorted raster plot (plotraster)
  """
  
  if len(dist)>0:
    k=lexsort([dist,classes])    
  else:
    k = classes.argsort();
  

  k=k.argsort()
  i=searchsorted(stim,spt);
  spt2=spt.copy();
  spt2=(spt2-stim[i-1]);
  spt3=spt2+stim[k[i-1]];
  spt3.sort()
  return spt3;

def kmedoids(nsteps,nclasses,s):
    """Looks for clusters in data with K-medoids
       Input: nsteps - maximal number of steps
              nclasses - number of classes
              s - observation matrix
     """
    
    (ntrials,ndim)=s.shape;
    medoids=random.randint(0,ntrials,nclasses)
    #dist=CalcDistance(s)
    member=ones(ntrials);
    for step in range(nsteps):
      member_new=s[:,medoids].argmin(1)
      if (member==member_new).all():
        break
      member=member_new.copy(); #Check if need to copy
      for i in range(nclasses):
        medoids[i]=sum(s[:,member==i],axis=1).argmin()
   
    w=0 
    for i in range(nclasses):
      Nk=sum(member==i)
      w+=sum(sum(s[:,member==i][member==i,:]))
   
    dist=s[:,medoids].min(1) 
    #return member,dist,w
    return member,dist
    
def kmeans(nsteps, nclasses, s):
    """Classifies vector by K-means algorithm
       Input: nsteps - maximal number of steps
              nclasses - number of classes
              s - observation matrix
       Output: 1-D array with class numbers
               D - clusters strength
    """
    (ntrials,ndim)=s.shape;
    centroids=rand(ntrials,nclasses)*(s.max()-s.min())+s.min();
    nearest=ones(ntrials);
    for step in range(nsteps):
      d=[];
      for i in range(nclasses):
        d.append(sum((s-centroids[:,i])**2,1));
      d=array(d);
      nearest_new=d.argmin(0);
      if (nearest==nearest_new).all():
        break
      nearest=nearest_new.copy();
      for i in range(nclasses):
        vecs,=where(nearest==i);
        centroids[:,i]=s[:,vecs].mean(1);
      
    #strength=zeros(nclasses,dtype='double');

    #w=0.0
    #for i in range(nclasses):
    #  Nk=sum(nearest==i)
    #  inclass,=where(nearest==i);
    #  outclass,=where(nearest<>i);
    #  strength[i]=sum(d[i,outclass])/sum(d[i,inclass])*double(Nk)/double(ntrials-Nk)
    #  w+=sum(d[i,inclass])*Nk;

    #dist=d.min(0)

    #nearest=nearest.astype(int) 
    #new_order=strength.argsort()[::-1]
    #nearest=new_order.argsort()[nearest]
    #return nearest,dist,strength,w
    return nearest

def CalcDistance(a):
    return array([[(sum((i-j)**2)) for i in a] for j in a])

def findK(nclasses,s,nrep=5):
    maxsteps=100;
    
    W=[];
    for i in nclasses:
      Wi=inf;
      for j in range(nrep):
        (cl,D,temp)=kmeans(maxsteps,i,s);
        if len(unique(cl))==i: Wi=min(temp,Wi)
      W.append(Wi);
    
    return W

def reshaping(s):
    (n,bins)=histogram(s.flatten(1),bins=50,normed=True);
    ncum=cumsum(concatenate(([0],n)))/50;
    bins=concatenate((bins,[1]))
    m=s.mean();
    (plsq,out)=optimize.leastsq(residuals, (m,0.1), args=(ncum, bins));
    y2=sigmoid(bins,plsq[0],plsq[1])
    #plt=plot(bins,y2,bins,ncum)
    return sigmoid(s,plsq[0],plsq[1])
    
def calcsimilarity(trains,width):
    len=100;
    y=gausswin(100,width);
    g=array([convolve(y,x) for x in trains]);
    norm=(sum(g**2,1));
    norm[norm==0]=1;
    gnorm=g/sqrt(norm[:,newaxis]);
    s=array([[dot(x,y) for x in gnorm] for y in gnorm]);
    return s,g

def BinTrains(trains,win,tsamp=0.25):
    """Convert a list of spike trains into a binary sequence"""
    #bins=arange(win[0]-tsamp, win[1]+2*tsamp, tsamp)
    bins=arange(win[0], win[1], tsamp)
    trains = [histogram(spt,bins)[0][1:] for spt in trains]
    #trains = [histogram(spt,bins)[0][1:-1] for spt in trains]

    return bins[1:-1],array(trains).T

def TimesToTrains(times,Fs=2E5,win=None):
    """Convert a list of spike times (in ms) to a sequence of 0s and 1s"""
    
    if not win:
      #Calculate window here
      win=[inf,-inf]
      for tr in times:
        win[0]=min(tr.min(),win[0])
        win[1]=max(tr.max(),win[1])

    npts=ceil((win[1]-win[0])/1000.0*Fs)+1
    t=zeros((len(times),npts))
    for i,x in enumerate(times):
      t[i,fix((x-win[0])/1000.0*Fs).astype(int)]=1
    return t 

def MetricRossumOld(times,tc):
    """Implements van Rossum metrics using convolution"""
    fs=2.E5
    #t=arange(0,100.0,dtype=float)/fs*1000
    tmax=3*log(10)*tc
    support=ceil(tmax/1000.0*fs);
    t=arange(0,support,dtype=float)/fs*1000
    trains=TimesToTrains(times,fs)
    y=exp(-t/tc)
    #y=y/sum(y)*tc
    g=array([convolve(x,y) for x in trains])
    s=array([[sum((tr1-tr2)**2) for tr1 in g] for tr2 in g])
    return s/tc/fs*1000

def MetricRossum(times,tc,norm=False):
    """Implements van Rossum metrics (van Rossum, Neural Computation 13,751 
 (2003))"""
    def R(tr1,tr2):
       """r=int(f(t)*g(t))dt=
           =tc/2 * sum_i sum_j (exp(-|t_i-t_j|/tc))"""
       #r=(sum(exp(tr1/tc))-sum(exp(tr2/tc)))**2
       r=sum(exp(-abs(subtract.outer(tr1,tr2).ravel())/tc))
       return r

    s_xx=array([R(x,x) for x in times])
    ntrials=len(times)
    s_xy=zeros((ntrials,ntrials))
    for i in xrange(ntrials):
        for j in xrange(i,ntrials):
            s_xy[i,j]=R(times[i],times[j])
    s_xy=s_xy+triu(s_xy,1).T
    s=s_xx[newaxis,:]+s_xx[:,newaxis]-2*s_xy

    if norm:
        nspikes=array([len(x) for x in times])
        s=s*2./(nspikes[:,newaxis]+nspikes[newaxis,:])
        s[isnan(s)]=0.
    return s/2
    
def MetricSchreiber(trains,width):
    """Calculate the similarity between spike trains, using the measure
    introduced by S. Schreiber & T. Sejnowski [Neurocomputing 52-54 (2003)]
    """
    sig42 = 4.0*width**2
    def f(s1,s2):
        if len(s1)<0 or len(s2)<0: return 0 
        d=subtract.outer(s1,s2).ravel() #differences of spike times
        return sum( exp(-d**2/sig42) )  #\int gauss(x,s1)*gauss(x,s2) dx
    
    #s=array([[f(x,y)/sqrt(f(x,x)*f(y,y)) for x in trains] for y in trains]);
    s=array([[f(x,y) for x in trains] for y in trains]);
    norm_coef=sqrt(diag(s)[:,newaxis]);
    s/=norm_coef;
    s/=norm_coef.T;
    s[isnan(s)]=0;
    return s

def MetricInterval(trains,q,append=True,norm=False):
   """Calculate similarity between spike trains based on interval metric as
   described by JD Victor & KP Purpura [Network 8:127-164 (1997)]"""
   if append:
     int=[diff(concatenate(([0],sp))) for sp in trains]
   else: 
     int=[diff(sp) for sp in trains]

   def G(s1,s2):
      """G_{i,j}=min(G_{i-1,j}+1,G_{i,j-1}+1,G_{i-1,j-1}+M(e_i,e_j))""" 
      if len(s1)==0 and len(s2)==0: return 0
      if len(s1)==0 or len(s2)==0: return 1
      

      g1=G(s1[:-1],s2)
      g2=G(s1,s2[:-1])
      g3=G(s1[:-1],s2[:-1]) + q*abs(s1[-1]-s2[-1])
      
      return min(g1+1,g2+1,g3)

   def Giter(s1,s2):
      m=zeros((len(s1)+1,len(s2)+1))
      m[0,:]=range(len(s2)+1); m[:,0]=range(len(s1)+1)
      for j in xrange(1,len(s2)+1):
         for i in xrange(1,len(s1)+1):
            m[i,j]=min(m[i-1,j]+1,m[i,j-1]+1,m[i-1,j-1]+q*abs(s1[i-1]-s2[j-1]))
      return m[-1,-1]

   if norm:
     s=array([[Giter(x,y)/(len(x)+len(y)) for x in int] for y in int]);
     s[isnan(s)]=0
   else: 
     s=array([[Giter(x,y) for x in int] for y in int]);
   return s   

def MetricSpike(trains,q,norm=False):
   """Calculate similarity between spike trains based on D^spike metric as
   described by JD Victor & KP Purpura [Network 8:127-164 (1997)]"""

   def Giter(s1,s2):
      m=zeros((len(s1)+1,len(s2)+1))
      #m[0,1:]=1; m[1:,0]=1
      m[0,:]=range(len(s2)+1); m[:,0]=range(len(s1)+1)
      for j in xrange(1,len(s2)+1):
         for i in xrange(1,len(s1)+1):
            m[i,j]=min(m[i-1,j]+1,m[i,j-1]+1,m[i-1,j-1]+q*abs(s1[i-1]-s2[j-1]))
      return m[-1,-1]
  
   if norm:
     s=array([[Giter(x,y)/(len(x)+len(y)) for x in trains] for y in trains]);
     s[isnan(s)]=0
   else: 
     s=array([[Giter(x,y) for x in trains] for y in trains]);
   return s   

def plotPSTH(spt, stim, win=[0,30], bin=0.25, ax=None, 
             rate=False, **kwargs):
   """Plot peri-stimulus time histogram (PSTH)"""
   i=searchsorted(stim+win[0],spt)
   spt2=(spt-stim[i-1])
   bins=arange(win[0],win[1],bin)
   (psth,bins)=histogram(spt2,bins)
   if not ax:
       ax=gca()
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
    ax=gca()
    time,binned=BinTrains(trains,win,bin)
    psth=mean(binned,1)/bin*1000
    #bar(time,psth,bin,**kwargs)
    plot(time,psth,**kwargs)
    if rate:
        psth=psth*1./len(stim)/bin*1000.
        ax.set_ylabel('firing rate (Hz)')
    else:
        ax.set_ylabel('number of spikes')
    ax.set_xlabel('time (ms)')
    ax.set_ylabel('firing rate (Hz)')

def CalcTrainsPSTH(trains, win, bin=0.25):
    time,binned=BinTrains(trains,win,bin)
    psth=mean(binned,1)/bin*1000
    return time, psth

def CalcPSTH(spt,stim,win=[0,30],bin=0.25,ax=None,norm=False,**kwargs):
   """Calculate peri-stimulus time histogram (PSTH).
   Output:
       -- psth - spike counts
       -- bins - bins edges"""
   i=searchsorted(stim+win[0],spt)
   spt2=(spt-stim[i-1])
   bins=arange(win[0],win[1],bin)
   (psth,bins)=histogram(spt2,bins)
   #(psth,bins)=histogram(spt2,bins)
   if norm:
       psth=psth*1000./(len(stim)*bin)
   return psth[1:],bins[1:-1]

def plotPSTHBar(spt,stim,win=[0,30],bin=0.25,**kwargs):
   """Plot peri-stimulus time histogram (PSTH)"""
   i=searchsorted(stim+win[0],spt)
   spt2=(spt-stim[i-1])
   bins=arange(win[0],win[1],bin)
   (psth,bins)=histogram(spt2, bins)
   psth=psth*1./len(stim)/bin*1000.
   gca().bar(bins[:-1], psth, bin, **kwargs)
   xlabel('time (ms)')
   ylabel('firing rate (Hz)')

#### UNIT TESTS ######

class PatternClassifier(unittest.TestCase):
    def testCountClasses_correct(self):
        a = np.array(['a','b','a','a','b'])
        n, l = CountClasses(a)
        self.failUnless((n==np.array([3,2])).all())
        self.failUnless((l==np.array(['a','b'])).all())
    def testCountClasses_correct(self):
        a = np.array(['a','b','a','a','b'])
        n, l = CountClasses(a,['b'])
        self.failUnless(n==2)

class SpikeMetricTests(unittest.TestCase):
    def testMetricRossum1(self):
        #symmetry test
        y=random.randint(0,15,(10,3))
        width=1.
        d=MetricRossum(y,width)
        d_norm=MetricRossum(y,width,norm=True)
        self.failUnless((triu(d)==tril(d).T).all())
        self.failUnless((triu(d_norm)==tril(d_norm).T).all())
    def testMetricRossum2(self):
        #non-degeneracy
        y=random.randint(0,15,10)
        width=0.01
        d=MetricRossum([y,y],width)
        d_norm=MetricRossum([y,y],width,norm=True)
        self.failUnless((d==zeros(2)).all())
        self.failUnless((d_norm==zeros(2)).all())
    def testMetricRossum3(self):
        #Limit cases tau->0 and tau->Inf
        y=[array([2.,3.,10.]),array([2.,5.,8.,10.])]
        d_0=MetricRossum(y,0.000001) #Coincidence detector
        d_inf=MetricRossum(y,10000.)
        self.assertAlmostEqual(1.500,d_0[0,1])
        self.assertAlmostEqual(0.500,d_inf[0,1],3)
    def testMetricRossum4(self):
        y=[array([ 7.6,  9.8,  4.4,  9.3,  1. ,  5.5,  7.4,  6.8]),
           array([ 5.7,  4.4,  5.2,  7.6,  1.2,  5. ,  2.5,  0.2,  3.1, 8.1])]
        d1=MetricRossum(y,0.5)
        d2=MetricRossum(y,1.)
        d3=MetricRossum(y,2.)
        self.assertAlmostEqual(6.48610928,d1[0,1])
        self.assertAlmostEqual(6.6687246,d2[0,1])
        self.assertAlmostEqual(7.07166336,d3[0,1])
        
        d1_norm=MetricRossum(y,0.5,norm=True)
        d2_norm=MetricRossum(y,1.,norm=True)
        d3_norm=MetricRossum(y,2.,norm=True)
        self.assertAlmostEqual(0.72067881,d1_norm[0,1])
        self.assertAlmostEqual(0.7409694,d2_norm[0,1])
        self.assertAlmostEqual(0.78574037,d3_norm[0,1])
    def testMetricRossum5(self):
        y=[array([2.,3.,10.]),array([4.,5.,8.,12.])]
        d_0=MetricRossum(y,0.000001,norm=True)
        self.assertAlmostEqual(1.000,d_0[0,1])

if __name__ == '__main__':
    unittest.main()
