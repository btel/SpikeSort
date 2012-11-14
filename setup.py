#!/usr/bin/env python
#coding=utf-8

import setuptools

from numpy.distutils.core import setup, Extension
 
diptst_ext = Extension(name = 'spike_sort.stats._diptst',
                       sources = ['src/spike_sort/stats/diptst/diptst.f',
                                  'src/spike_sort/stats/diptst/diptst.pyf'])


setup(name='SpikeSort',
      version='0.1',
      description='Python Spike Sorting Package',
      author='Bartosz Telenczuk',
      author_email='bartosz.telenczuk@gmail.com',
      url='http://neuroscience.telenczuk.pl',

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


