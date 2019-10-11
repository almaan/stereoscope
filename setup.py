#!/usr/bin/env python3

from setuptools import setup
import os

setup(name='STereoSCope',
            version='0.1',
            description='Integration of ST and SC data',
            url='http://github.com/almaan/stsc',
            author='Alma Andersson',
            author_email='alma.andersson@scilifelab.se',
            license='MIT',
            packages=['STereoSCope'],
            install_requires=[
                            'torch',
                            'numpy',
                            'pandas',
                            'scipy',
                            'logging',
                            'argparse',
                            'matplotlib',
                            'sklearn',
                            'umap-learn',
                      ],
            zip_safe=False)



