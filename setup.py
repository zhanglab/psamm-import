#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='psamm-import',
    version='0.1',
    description='PSAMM model importers',
    maintainer='Jon Lund Steffensen',
    maintainer_email='jon_steffensen@uri.edu',
    url='https://github.com/zhanglab/psamm-import',

    packages=find_packages(),
    entry_points = {
        'console_scripts': [
            'psamm-import = psamm_import.importer:main'
        ]
    },

    install_requires=[
        'PyYAML>=3.11,<4.0',
        'xlrd',
        'psamm'
    ])
