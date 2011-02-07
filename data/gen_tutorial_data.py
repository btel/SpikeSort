#!/usr/bin/env python
#coding=utf-8

from spike_sort.io.filters import PyTablesFilter, BakerlabFilter

in_dataset = "/Gollum/s5gollum01/el3"
out_dataset = "/SubjectA/session01/el1/raw"

in_filter = BakerlabFilter("gollum.inf")
out_filter = PyTablesFilter("tutorial.h5")

sp = in_filter.read_sp(in_dataset)
out_filter.write_sp(sp, out_dataset)

in_filter.close()
out_filter.close()
