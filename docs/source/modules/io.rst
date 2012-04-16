File I/O (:mod:`spike_sort.io`)
===============================

.. currentmodule:: spike_sort.io

Functions for reading and writing datafiles.

.. _io_filters:

Read/Write Filters (:mod:`spike_sort.io.filters`)
-------------------------------------------------

Filters are basic backends for read/write operations. They offer following 
methods:

 * :py:func:`read_spt` -- read event times (such as spike times)
 * :py:func:`read_sp` -- read raw spike waveforms 
 * :py:func:`write_spt` -- write spike times
 * :py:func:`write_sp` -- write raw spike waveforms

The `read_*` methods take usually one argument (`datapath`), but it is
not required. The `write_*` methods take `datapath` and the data to be
written. 

If you want to read/write you custom data format, it is enough that
you implement a class with these functions.

The following filters are implemented:

.. autosummary::

   filters.BakerlabFilter
   filters.PyTablesFilter


Export tools (:mod:`spike_sort.io.export`)
--------------------------------------------

These tolls take one of the `io.filters` as an argument and export
data to the file using `write_spt` or `write_sp` methods.

.. autosummary::

   export.export_cells

Reference
---------

.. automodule:: spike_sort.io.filters
   :members: 

.. automodule:: spike_sort.io.export
   :members: 

