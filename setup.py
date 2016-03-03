#!/usr/bin/env python
#coding=utf-8

import setuptools

from numpy.distutils.core import setup, Extension

ext_modules = []

import os
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:
    diptst_ext = Extension(name = 'spike_sort.stats._diptst',
                        sources = ['src/spike_sort/stats/diptst/diptst.f',
                                  'src/spike_sort/stats/diptst/diptst.pyf'])
    ext_modules.append(diptst_ext)


setup(name='SpikeSort',
      version='0.13',
      description='Python Spike Sorting Package',
      long_description="""SpikeSort is a flexible spike sorting framework 
      implemented entirely  in Python based on widely-used packages such as numpy,
      PyTables and matplotlib. It features manual and automatic clustering, many 
      data formats and it is memory-efficient.""",
      author='Bartosz Telenczuk and Dmytro Bielievtsov',
      author_email='bartosz.telenczuk@gmail.com',
      url='http://spikesort.org',
      ext_modules = ext_modules,
      classifiers = [
          "Development Status :: 4 - Beta",
          "Environment :: Console",
          "Environment :: X11 Applications",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: BSD License",
          "Operating System :: OS Independent",
          "Programming Language :: Python :: 2.7"
          ],
      packages=['spike_sort', 
                'spike_sort.core', 
                'spike_sort.stats',
                'spike_sort.ui',
                'spike_sort.io', 
                'spike_beans',
                'spike_analysis',
                ],
      package_dir = {"": "src"},
      install_requires=[
          'matplotlib',
          'tables',
          'numpy >= 1.4.1',
          'scipy',
          'PyWavelets'
        ]
)


