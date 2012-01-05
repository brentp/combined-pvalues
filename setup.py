#!/usr/bin/env python

from distutils.core import setup

setup(name='cpv',
      version='0.1',
      description='combine p-values',
      author='Brent Pedersen',
      author_email='bpederse@gmail.com',
      license='MIT',
      url='https://github.com/brentp/combined-pvalues',
      packages=['cpv', 'cpv.tests'],
      scripts=['cpv/comb-p'],
      long_description=open('README.rst').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
 )
