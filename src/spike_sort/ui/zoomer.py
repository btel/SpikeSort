class Zoomer(object):
    '''Allows to zoom subplots'''
    def __init__(self, plt, fig):
        self.axlist=[]
        
        self.zoomed_state={'geometry':(1,1,1),
                           'xlabel_visible':True,
                           'ylabel_visible':True}
        
        self.old_state={'geometry':[1,1,1],
                        'xlabel_visible':True,
                        'ylabel_visible':True}
        
        plt.connect('key_press_event', self.zoom)
        self.fig=fig
    
    def zoom(self, event):
        axis=event.inaxes
        
        if axis==None or event.key!='z':
            return
        
        if axis.get_geometry() == self.zoomed_state['geometry']:
            zoomed=True
        else:
            zoomed=False
        
        if zoomed==False:
            #saving previous state
            self.old_state['geometry']=axis.get_geometry()
            self.old_state['xlabel_visible']=axis.xaxis.label.get_visible()
            self.old_state['ylabel_visible']=axis.yaxis.label.get_visible()
            
            #removing old axes
            self.axlist=list(self.fig.get_axes())
            for ax in self.axlist:
                self.fig.delaxes(ax)
            
            #modifying state (zooming)
            axis.change_geometry(*self.zoomed_state['geometry'])
            axis.xaxis.label.set_visible(self.zoomed_state['xlabel_visible'])
            axis.yaxis.label.set_visible(self.zoomed_state['ylabel_visible'])
            self.fig.add_axes(axis)
            self.fig.show()
            
        else:
            #removing current axes
            self.fig.delaxes(axis)
            
            #bringing the old state back
            axis.change_geometry(*self.old_state['geometry'])
            axis.xaxis.label.set_visible(self.old_state['xlabel_visible'])
            axis.yaxis.label.set_visible(self.old_state['ylabel_visible'])
            
            #adding old axes
            for ax in self.axlist:
                self.fig.add_axes(ax)
            
            self.fig.show()
