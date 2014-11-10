#!/usr/bin/env python
import ez_setup
ez_setup.use_setuptools()


from setuptools import setup
import cpv
version = cpv.__version__

setup(name='cpv',
      version=version,
      description='combine p-values',
      author='Brent Pedersen',
      author_email='bpederse@gmail.com',
      license='MIT',
      url='https://github.com/brentp/combined-pvalues',
      packages=['cpv', 'cpv.tests'],
      install_requires=['scipy', 'numpy', 'toolshed', 'interlap'],
      scripts=['cpv/comb-p'],
      long_description=open('README.rst').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
 )
