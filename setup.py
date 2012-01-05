#!/usr/bin/env python

from distutils.core import setup

setup(name='combp',
      version='0.1',
      description='combine p-values',
      author='Brent Pedersen',
      author_email='bpederse@gmail.com',
      url='https://github.com/brentp/combined-pvalues',
      packages=['cpv', 'cpv.tests'],
      scripts=['cpv/comb-p']
     )
