#!/usr/bin/env python3

from setuptools import setup
import os

setup(name='STereoSCope',
            version='0.2',
            description='Integration of ST and SC data',
            url='http://github.com/almaan/stsc',
            author='Alma Andersson',
            author_email='alma.andersson@scilifelab.se',
            license='MIT',
            packages=['stsc'],
            install_requires=[
                            'torch>=1.1.0',
                            'numba',
                            'numpy>=1.14.0',
                            'pandas>=0.25.0',
                            'logging',
                            'argparse',
                            'matplotlib>=3.1.0',
                            'sklearn',
                            'umap-learn>=0.3.10',
                            'Pillow',
                      ],
            entry_points={'console_scripts': ['stereoscope = stsc.__main__:main',
                                             ]
                         },
            zip_safe=False)

