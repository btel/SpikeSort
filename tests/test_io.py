import spike_sort as ss

from nose.tools import ok_, eq_, raises
import tables
import numpy as np
import os
import spike_sort.io.hdf5
import spike_sort.io.bakerlab
import filecmp
import json

class TestHDF:
    def setUp(self):
        self.data = np.random.randint(1000,size=(100,4))
        self.spt = np.random.randint(0,100, (10,))/200.
        self.el_node = '/Subject/Session/Electrode'
        self.fname = 'test.h5'
        self.cell_node = self.el_node+'/cell'
        self.h5f =  tables.openFile(self.fname,'a')
        
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
        
class TestBakerlab:
    def setup(self):
        file_descr = {"fspike":"{ses_id}{el_id}.sp",
                      "dirname":".",
                      "FS":5.E3}
        self.el_node = '/Test/s32test01/el1'
        self.data = np.random.randint(-1000, 1000, (100,))/200.
        self.spt = np.random.randint(0,100, (10,))
        self.conf_file = 'test.conf'
        self.fname = "32test011.sp"
        
        with open(self.conf_file, 'w') as fp:
             json.dump(file_descr, fp)
          
        (self.data*200).astype(np.int16).tofile(self.fname)
        
    def tearDown(self):
        os.unlink(self.conf_file)
        os.unlink(self.fname)
        
    def test_write_spt(self):
        pass
    def test_read_spt(self):
        
        pass
    def test_write_sp(self):
        el_node_tmp = '/Test/s32test01/el2'
        spike_sort.io.bakerlab.write_sp({"data":self.data}, self.conf_file,
                                       el_node_tmp)
        files_eq = filecmp.cmp("32test011.sp","32test012.sp", shallow=0)
        os.unlink("32test012.sp")
        ok_(files_eq)
        
    def test_read_sp(self):
        sp = spike_sort.io.bakerlab.read_sp(self.conf_file,
                                            self.el_node)
        read_data = sp['data']
        ok_((np.abs(read_data-self.data)<=1/200.).all())
        