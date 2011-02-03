'''
Created on Feb 2, 2011

@author: bartosz
'''

from spike_sort import cluster

def export_cells(io_filter, node_templ, spike_times, overwrite=False):
    """
    Export discriminated spike times of all cells to a file.
    
    :arguments:
     * io_filter -- read/write filter object (see :py:mod:`spike_sort.io.filters`)
     * node_templ -- template (format string) for the dataset identifiers
     * spt_dict -- mapping object with keys: :py:attr:`data` (spike times) and 
       :py:attr:`cell_id` (cluster/cell id number)
    """

    #spike_times = cluster.split_cells(spt_dict, spt_dict['cell_id'])
    
    for cell_id, spt_cell in spike_times.items():
        dataset = node_templ.format(cell_id=cell_id)
        io_filter.write_spt(spt_cell, dataset, overwrite=overwrite)