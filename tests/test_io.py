import spike_sort as ss

from nose.tools import ok_, eq_, raises
import tables
import numpy as np
import os
import spike_sort.io.hdf5
import filecmp

class TestHDF:
    def setUp(self):
        self.fname = 'test.h5'
        self.el_node = '/Subject/Session/Electrode'
        self.cell_node = self.el_node+'/cell'
        self.h5f =  tables.openFile(self.fname,'a')
        self.data = np.random.rand(100,4)
        self.spt = np.random.randint(0,100, (10,))
        self.spt.sort()
        atom = tables.Atom.from_dtype(self.data.dtype)
        shape = self.data.shape
        filter = tables.Filters(complevel=0, complib='zlib')
        new_array = self.h5f.createCArray(self.el_node, "raw", atom, shape, 
                                     filters=filter, 
                                     createparents=True)
        self.sampfreq=5.E3
        new_array.attrs['sampfreq']=self.sampfreq
        new_array[:] = self.data
        spt_array = self.h5f.createArray(self.el_node, "cell", self.spt, 
                                         title="",
                                         createparents="True")
        self.h5f.close()
        
    def tearDown(self):
        spike_sort.io.hdf5.close_all()
        os.unlink(self.fname)
    
    def test_write(self):
        sp_dict = {'data':self.data,'FS':self.sampfreq}
        spt_dict = {'data':self.spt}
        spike_sort.io.hdf5.write_sp(sp_dict, "test2.h5", self.el_node+"/raw")
        spike_sort.io.hdf5.write_spt(spt_dict, "test2.h5", self.cell_node)
        spike_sort.io.hdf5.close_all()
        exit_code = os.system('h5diff ' + self.fname + ' test2.h5')
        os.unlink("test2.h5")
        ok_(exit_code==0)
        
    def test_read_sp(self):
        sp = spike_sort.io.hdf5.read_sp(self.fname, self.el_node)
        ok_((sp['data'][:]==self.data).all())
    
    def test_read_spt(self):
        spt = spike_sort.io.hdf5.read_spt(self.fname, self.cell_node)
        ok_((spt['data']==self.spt).all())