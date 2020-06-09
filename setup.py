#!/usr/bin/env python3

from setuptools import setup
import os
import sys

setup(name='STereoSCope',
            version='0.2',
            description='Integration of ST and SC data',
            url='http://github.com/almaan/stereoscope',
            author='Alma Andersson',
            author_email='alma.andersson@scilifelab.se',
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
                            'Pillow',
                      ],
            entry_points={'console_scripts': ['stereoscope = stsc.__main__:main',
                                             ]
                         },
            zip_safe=False)

