def export_cells(io_filter, node_templ, spike_times, overwrite=False):
    """Export discriminated spike times of all cells to a file.
    
    Parameters
    ----------
    io_filter : object,
        read/write filter object (see :py:mod:`spike_sort.io.filters`)
    node_templ : string
        string identifing the dataset name. It will be passed to
        IOFilters.write_spt method. It can contain the
        `{cell_id}` placeholder that will be substituted by cell
        identifier. 
    spt_dict : dict
        dictionary in which keys are the cell IDs and values are spike
        times structures
    """
    for cell_id, spt_cell in spike_times.items():
        dataset = node_templ.format(cell_id=cell_id)
        io_filter.write_spt(spt_cell, dataset, overwrite=overwrite)
