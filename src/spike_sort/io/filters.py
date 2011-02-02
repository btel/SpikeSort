'''
Created on Feb 2, 2011

@author: bartosz
'''
import os
import numpy as np
import json
import re
from tempfile import mkdtemp
import tables

class BakerlabFilter:
    
    def __init__(self, conf_file):
        self._regexp="^/(?P<subject>[a-zA-z]+)/s(?P<ses_id>.+)/el(?P<el_id>[0-9]+)(/cell)?(?P<cell_id>[0-9]+)?$"
        self.conf_file = conf_file
        with file(conf_file) as fid:
            self.conf_dict = json.load(fid)
            

    
    def read_sp(self, dataset, memmap=None):
        """Reads raw spike waveform from file in bakerlab format
        
        :arguments:
         - conf_file -- configuration file
         - dataset -- dataset path
         - memmap -- use memory mapped arrays to save some memory
        """
       
        conf_dict = self.conf_dict
        m = re.match(self._regexp, dataset)
        rec_dict = m.groupdict()
        n_contacts = conf_dict['n_contacts']
        f_spike = conf_dict['fspike']
        
        dirname = conf_dict['dirname']
        rec_dict['contact_id']=1
        
        full_path = os.path.join(dirname, f_spike)
        fname = full_path.format(**rec_dict)
        sp = np.fromfile(fname, dtype=np.int16)
    
        
        if memmap=="numpy":
            #create temporary memmory mapped array
            filename = os.path.join(mkdtemp(), 'newfile.dat')
            fp = np.memmap(filename, dtype='float', mode='w+', 
                           shape=(len(sp),n_contacts))
        elif memmap=="tables":
            atom = tables.Atom.from_dtype(sp.dtype)
            shape = (len(sp), n_contacts)
            filters = tables.Filters(complevel=3, complib='blosc')
            filename = os.path.join(mkdtemp(), 'newfile.dat')
            h5f = tables.openFile(filename,'w')
            fp = h5f.createCArray('/', "test", atom, shape, filters=filters)
        else:
            fp = np.empty((len(sp), n_contacts), dtype=sp.dtype)
    
        fp[:,0]=sp
        for i in range(1,n_contacts):
            rec_dict['contact_id']=i+1
            fname = full_path.format(**rec_dict)
            sp = np.fromfile(fname,dtype=np.int16)
            fp[:,i]=sp
        del sp
        return {'data':fp, "FS":conf_dict['FS'], "n_contacts":n_contacts} 
    
    def write_sp(self, sp_dict, dataset):
        """Write raw spike waveform to a file in bakerlab format
        
        :arguments:
         - sp_dict -- spike waveform dict
         - conf_file -- configuration file
         - dataset -- dataset path
        """
        sp = sp_dict['data']
        conf_dict = self.conf_dict
    
        m = re.match(self._regexp, dataset)
        n_contacts = sp.shape[1]
        rec_dict = m.groupdict()
        
        for i in range(n_contacts):
            rec_dict['contact_id']=i+1
            fname = conf_dict['fspike'].format(**rec_dict)
            full_path = os.path.join(conf_dict['dirname'], fname)
            sp_int = (sp[:,i]).astype(np.int16)
            sp_int.tofile(full_path)
    
    def _match_dataset(self, dataset):
        m = re.match(self._regexp, dataset)
        if not m:
            StandardError("dataset id could not be parsed")
        return m.groupdict()
    
    def read_spt(self,  dataset):
        """Returns spike times in miliseconds:
        
        :arguments:
         * dir_name : directory names with the data
         * dataset : dataset path
        """
        conf_dict = self.conf_dict
        rec = self._match_dataset(dataset)
        
        fname = conf_dict['fspt'].format(**rec)
        spt = np.fromfile(fname, dtype=np.int32)
        return {"data": spt/200.}
    
    
    def write_spt(self, spt_dict, dataset):
        """Returns spike times in miliseconds:
        
        :Arguments:
         * dir_name : directory names with the data
         * dataset : dataset name
        """
        
        conf_dict = self.conf_dict
        rec = self._match_dataset(dataset)
        spt = spt_dict['data']
   
        fname = conf_dict['fspt'].format(**rec)
        export_spt = (spt*200).astype(np.int32)
        export_spt.tofile(fname)
        
class PyTablesFilter:
    """
    HDF5 is a hierarchical datafile -- data is organised in a tree. The
    standard layout is::
       
        /{SubjectName}/
        /{SubjectName}/{SessionName}/{ElectrodeID}/
        /{SubjectName}/{SessionName}/{ElectrodeID}/stim: stimulus time
        /{SubjectName}/{SessionName}/{ElectrodeID}/raw: spike waveforms
        /{SubjectName}/{SessionName}/{ElectrodeID}/{CellID}: spike waveforms
        /{SubjectName}/{SessionName}/{ElectrodeID}/{CellID}/spt: spike
        times
    
    where curly brackets `{}` denote a group.
    
    This layout may be adjusted by changing paths
    
    TODO: implement custom layouts
    """
    
    def __init__(self, fname, mode='a'):
        self.h5file = tables.openFile(fname, mode) 
    
    @staticmethod
    def _get_attrs(node):
        PYTABLES_ATTRS = ["VERSION", "TITLE", "FLAVOR", "CLASS"]
        extra_attrs=dict([(name, node.getAttr(name)) 
            for name in node.attrs._v_attrnames
            if name not in PYTABLES_ATTRS])
    
        return extra_attrs
    
    def _open_file(self, fname, mode='r'):
        if type(fname) is tables.File:
            h5f = fname
        else:
            try:
                h5f = tables.openFile(fname, mode)
            except ValueError:
                fname.close()
                _open_files.remove(fname)
                h5f = tables.openFile(fname, mode)
        if not h5f in _open_files:
            _open_files.append(h5f)
        return h5f
    
    def close_all(self):
        global _open_files
        for f in _open_files:
            f.close()
        _open_files=[]
    
    def read_sp(self, dataset):
        """Read continous waveforms (EEG, LFG, spike waveform.
        
        :arguments:
            * dataset -- (string) path pointing to cell node
        """
        
        h5f = self.h5file
        
        electrode = "/".join(dataset.split('/')[:4])
    
        electrode_node = h5f.getNode(electrode)
    
        sp_raw = electrode_node.raw
        FS = electrode_node.raw.attrs['sampfreq']
    
        try:
            n_contacts = sp_raw.shape[1]
        except IndexError:
            n_contacts = 1
        
        return {"data": sp_raw, "FS": FS, "n_contacts":n_contacts} 
    
    def read_spt(self,  dataset):
        """Read event times (such as spike or stimulus times).
       
        :arguments:
            * dataset -- (string) with path pointing to cell node
        """
    
        h5f = self.h5file
    
        cell_node = h5f.getNode(dataset)
    
        spt = cell_node.read()[:].copy()
        
        ret_dict = {"data" :spt}
        extra_attrs = self._get_attrs(cell_node)
        ret_dict.update(extra_attrs)
       
        cell_node.flush()
        #h5f.close()
        
        return ret_dict 
    
    def write_spt(self, spt_dict, dataset,overwrite=False):
        
        h5f = self.h5file
        
        spt = spt_dict['data']
        
        parts = dataset.split('/')
        group = '/'.join(parts[:-1])
        node_name = parts[-1]
    
        if overwrite:
            try:
                h5f.removeNode(group, node_name)
            except tables.exceptions.NodeError:
                pass
    
        arr_node = h5f.createArray(group, node_name, spt, 
                title="" , createparents=True)
    
        attrs = spt_dict.copy()
        del attrs['data']
    
        for k, v in attrs.items():
            arr_node.setAttr(k, v)
    
        #arr_node._v_attrs.__dict__.update(attrs)
        
        #arr_node.flush()
        #h5f.close()
        
    def write_sp(self, sp_dict, dataset,overwrite=False):
        
        h5f = self.h5file
        
        sp = sp_dict['data']
        
        parts = dataset.split('/')
        group = '/'.join(parts[:-1])
        node_name = parts[-1]
    
        if overwrite:
            try:
                h5f.removeNode(group, node_name)
            except tables.exceptions.NodeError:
                pass
        
        atom = tables.Atom.from_dtype(sp.dtype)
        shape = sp.shape
        filters = tables.Filters(complevel=0, complib='zlib')
        arr_node = h5f.createCArray(group, node_name, atom, shape, 
                                         filters=filters, 
                                         createparents=True)
        arr_node[:] = sp
        
        arr_node.attrs['sampfreq']=sp_dict['FS']
        #attrs = sp_dict.copy()
        #del attrs['data']
        #del attrs['FS']
    
        #for k, v in attrs.items():
        #    arr_node.setAttr(k, v)
    
        #arr_node._v_attrs.__dict__.update(attrs)
        
        #arr_node.flush()
        #h5f.close()
    def close(self):
        if self.h5file:
            self.h5file.close()
            self.h5file = None

