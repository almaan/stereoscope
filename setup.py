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
                            'torch',
                            'numpy',
                            'pandas',
                            'logging',
                            'argparse',
                            'matplotlib',
                            'sklearn',
                            'umap-learn',
                            'Pillow',
                      ],
            entry_points={'console_scripts': ['stereoscope = stsc.__main__:main',
                                             ]
                         },
            zip_safe=False)

