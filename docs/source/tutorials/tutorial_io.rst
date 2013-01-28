.. _io_tutorial:

Reading custom data formats
===========================

.. testsetup::
   
.. testcleanup::

SpikeSort provides a flexible interface with various data sources via
so call input-output `Filters`. The code available for download already
contains a few filters that allow for reading (and in some cases writing)
of several data formats, including:

* raw binary data,
* HDF5 files,
* proprietary data formats (such as Axon `.abf`) supported by 
  `Neo <http://packages.python.org/neo/>`_ library.

In this tutorial we show you how to use available filters to read and plot data
from `.abf` files and next we will convience you how easy it is to add support 
for custom formats by defining your own `Filter`.

   

