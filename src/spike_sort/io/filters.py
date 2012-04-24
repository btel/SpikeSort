'''
Filters are basic backends for read/write operations. They offer following
methods:

 * read_spt -- read event times (such as spike times)
 * write_spt -- write spike times
 * read_sp -- read raw spike waveforms
 * write_sp -- write raw spike waveforms

'''

import os
import numpy as np
import json
import re
from tempfile import mkdtemp
import tables
import os.path


class BakerlabFilter(object):
    """Filter for custom binary data structure.

    The binary data consists of independent files for each contact in
    each electrode written as 16-bit ints.

    The paths to the datafiles are defined in :file:`.inf` file that is
    JSON-compatible and contains at least the following attributes:

      fspike : str
          path to raw recordings relative to `dirname`
      cell : str
          path to spike times (with resolution 20 us) relative to
          `dirname`
      n_contacts : int
          number of contacts per electrode
      dirname : str
          path to the data
      FS : int
          spike sampling frequency

    Each of the paths can include any of the following Python formatting
    placeholders:

      * `{subject}` -- subeject name
      * `{cell_id}` -- cell id
      * `{ses_id}` -- session id
      * `{el_id}` -- electrode id

    These will be substituted by data extracted from `datapath`
    paramaeters of the reading and writing methods.

    Parameters
    ----------

    conf_file : str
        path to the configuration file

    """
    def __init__(self, conf_file):
        """constructor"""
        self._regexp = "^/(?P<subject>[a-zA-z]+)/s(?P<ses_id>.+)/el(?P<el_id>[0-9]+)/?(?P<type>[a-zA-Z]+)?(?P<cell_id>[0-9]+)?$"
        self.conf_file = conf_file
        self._tempfiles = []
        with open(conf_file) as fid:
            self.conf_dict = json.load(fid)
            self.chunksize = 10E6  # number of elements in a chunk

    def read_sp(self, dataset, memmap=None):
        """Reads raw spike waveform from file in bakerlab format

        Parameters
        ----------
        dataset : str
            dataset path (in format
            /{subject}/session{ses_id}/el{el_id})
        memmap : {'numpy', 'tables', None}, optional
            if True use memory mapped arrays to save some memory
            (defaults to no memmory-mapping)
        """
        conf_dict = self.conf_dict
        m = re.match(self._regexp, dataset)
        rec_dict = m.groupdict()
        n_contacts = conf_dict['n_contacts']
        f_spike = conf_dict['fspike']

        dirname = conf_dict['dirname'].format(**os.environ)
        rec_dict['contact_id'] = 1

        full_path = os.path.join(dirname, f_spike)
        fname = full_path.format(**rec_dict)
        npts = os.path.getsize(fname) / 2
        # sp = np.memmap(fname, dtype=np.int16, mode='r+')
        dtype = 'int16'
        shape = (n_contacts, npts)

        if memmap == "numpy":
            # create temporary memory mapped array
            filename = os.path.join(mkdtemp(), 'newfile.dat')
            fp = np.memmap(filename, dtype=np.int16, mode='w+', shape=shape)
            self._tempfiles.append(fp)
        elif memmap == "tables":
            atom = tables.Atom.from_dtype(np.dtype(dtype))
            filters = tables.Filters(complevel=0, complib='blosc')
            filename = os.path.join(mkdtemp(), 'newfile.dat')
            h5f = tables.openFile(filename, 'w')
            self._tempfiles.append(h5f)
            fp = h5f.createCArray('/', "test", atom, shape, filters=filters)
        else:
            fp = np.empty(shape, dtype=np.int16)

        sz = np.min([self.chunksize, npts])
        n_chunks = int(np.ceil(npts / sz))
        for i in range(n_contacts):
            rec_dict['contact_id'] = i + 1
            fname = full_path.format(**rec_dict)
            sp = np.memmap(fname, dtype=np.int16, mode='r')
            
            # read and copy data by chunks
            for j in range(n_chunks):
                stop = np.min((npts, (j + 1) * sz))
                fp[i, j * sz:stop] = sp[j * sz:stop]
            # fp[:,i]=sp[:]
            del sp
        return {'data': fp, "FS": conf_dict['FS'], "n_contacts": n_contacts}

    def write_sp(self, sp_dict, dataset):
        """Write raw spike waveform to a file in bakerlab format

        Parameters
        ----------
        sp_dict : dict
            spike waveform dict
        dataset : str
            dataset path

        See Also
        --------
        read_sp
        """
        sp = sp_dict['data']
        conf_dict = self.conf_dict

        m = re.match(self._regexp, dataset)
        n_contacts = sp.shape[0]
        rec_dict = m.groupdict()
        dirname = conf_dict['dirname'].format(**os.environ)

        for i in range(n_contacts):
            rec_dict['contact_id'] = i + 1
            fname = conf_dict['fspike'].format(**rec_dict)
            full_path = os.path.join(dirname, fname)
            sp_int = (sp[i, :]).astype(np.int16)
            sp_int.tofile(full_path)

    def _match_dataset(self, dataset):
        m = re.match(self._regexp, dataset)
        if not m:
            raise StandardError("dataset id could not be parsed")
        return m.groupdict()

    def read_spt(self,  dataset):
        """Returns spike times in miliseconds:

        Parameters
        ----------
        dataset : str
            dataset path in format
            /{subject}/session{ses_id}/el{el_id}/cell{cell_id}
        """
        conf_dict = self.conf_dict
        rec = self._match_dataset(dataset)

        dirname = conf_dict['dirname'].format(**os.environ)
        fspt = conf_dict[rec['type']]
        full_path = os.path.join(dirname, fspt)
        fname = full_path.format(**rec)

        spt = np.fromfile(fname, dtype=np.int32)
        return {"data": spt / 200.0}

    def write_spt(self, spt_dict, dataset, overwrite=False):
        """Returns spike times in miliseconds.

        Parameters
        ----------
        dataset : str
            dataset path

        See Also
        --------
        read_spt
        """

        conf_dict = self.conf_dict
        rec = self._match_dataset(dataset)
        spt = spt_dict['data']
        dirname = conf_dict['dirname'].format(**os.environ)

        fspt = conf_dict[rec['type']]
        full_path = os.path.join(dirname, fspt)
        fname = full_path.format(**rec)
        if os.path.exists(fname) and not overwrite:
            raise IOError("file {0} already exists".format(fname))
        export_spt = (spt * 200).astype(np.int32)
        export_spt.tofile(fname)

        if 'metadata' in spt_dict:
            log_fname = fname[:-3] + 'log'
            if os.path.exists(log_fname) and not overwrite:
                raise IOError("file {0} already exists".format(log_fname))
            logfile = open(log_fname, 'w')
            json.dump(spt_dict['metadata'], logfile)
            logfile.close()

    def close(self):
        for f in self._tempfiles:
            fname = f.filename
            del f
            os.unlink(fname)
        self._tempfiles = []


class PyTablesFilter(object):
    """
    Read/Write HDF5

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
    """

    def __init__(self, fname, mode='a'):
        self.h5file = tables.openFile(fname, mode)

    @staticmethod
    def _get_attrs(node):
        PYTABLES_ATTRS = ["VERSION", "TITLE", "FLAVOR", "CLASS"]
        extra_attrs = dict([(name, node.getAttr(name))
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
        _open_files = []

    def read_sp(self, dataset):
        """Read continous waveforms (EEG, LFG, spike waveform)

        Parameters
        ----------
        dataset : str
            path pointing to cell node
        """

        h5f = self.h5file

        electrode = '/'.join(dataset.split('/')[:4])

        electrode_node = h5f.getNode(electrode)

        sp_raw = electrode_node.raw
        FS = electrode_node.raw.attrs['sampfreq']

        try:
            n_contacts, _ = sp_raw.shape
        except ValueError:
            n_contacts = 1
            sp_raw = sp_raw.read()[np.newaxis, :]

        return {"data": sp_raw, "FS": FS, "n_contacts": n_contacts}

    def read_spt(self,  dataset):
        """Read event times (such as spike or stimulus times).

        Parameters
        ----------
        dataset : str
            path pointing to cell node
        """

        h5f = self.h5file

        cell_node = h5f.getNode(dataset)

        spt = cell_node.read()[:].copy()

        ret_dict = {"data": spt}
        extra_attrs = self._get_attrs(cell_node)
        ret_dict.update(extra_attrs)

        cell_node.flush()
        #h5f.close()

        return ret_dict

    def write_spt(self, spt_dict, dataset, overwrite=False):
        """Write spike times"""
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
                title="", createparents=True)

        attrs = spt_dict.copy()
        del attrs['data']

        for k, v in attrs.items():
            arr_node.setAttr(k, v)

    def write_sp(self, sp_dict, dataset, overwrite=False):
        """Write signal"""
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

        arr_node.attrs['sampfreq'] = sp_dict['FS']

    def close(self):
        if self.h5file:
            self.h5file.close()
            self.h5file = None
