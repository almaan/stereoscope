#!/usr/bin/env python3

from setuptools import setup
import os
import sys

setup(name='stereoscope',
            version='0.3',
            description='Integration of ST and SC data',
            author='Alma Andersson',
            author_email='alma.andersson@scilifelab.se',
            url='http://github.com/almaan/stereoscope',
            download_url='https://github.com/almaan/stereoscope/archive/v_03.tar.gz',
            license='MIT',
            packages=['stsc'],
            python_requires='>3.0.0',
            install_requires=[
                            'torch>=1.1.0',
                            'numba>=0.49.0',
                            'numpy>=1.14.0',
                            'pandas>=0.25.0',
                            'matplotlib>=3.1.0',
                            'scikit-learn>=0.20.0',
                            'umap-learn>=0.4.1',
                            'anndata',
                            'scipy',
                            'Pillow',
                      ],
            entry_points={'console_scripts': ['stereoscope = stsc.__main__:main',
                                             ]
                         },
            zip_safe=False)

