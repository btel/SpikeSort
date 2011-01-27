#!/usr/bin/env python
#coding=utf-8

from setuptools import setup

setup(name='SpikeSort',
      version='0.1',
      description='Python Spike Sorting Package',
      author='Bartosz Telenczuk',
      author_email='bartosz.telenczuk@gmail.com',
      url='http://neuroscience.telenczuk.pl',
      packages=['spike_sort', 'spike_sort.core', 'spike_sort.ui',
                'spike_sort.io'],
      package_dir = {"": "src"},
      test_suite='nose.collector'
     )

