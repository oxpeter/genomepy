#!/usr/bin/env python

from distutils.core import setup

setup(name='Genomepy',
      version='1.0',
      description='Python Distribution Utilities',
      author='Peter Oxley',
      author_email='oxpeter+git@gmail.com',
      url='https://github.com/oxpeter/genomepy',
      packages=['genomepy'],
      package_dir={'genomepy': 'genomepy'},
      package_data={'genomepy': ['data/*.cfg'], 'genomepy': ['data/README']},
      requires=['argparse','scipy','numpy','matplotlib','Bio','qvalue', 'progressbar',
                'operator', 'pylab', 'mpl_toolkits', 'statsmodels', 'pysam'
                ]
     )