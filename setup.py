#!/usr/bin/env python
# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

from setuptools import setup, find_packages

setup(
    name='psamm-import',
    version='0.3',
    description='PSAMM model importers',
    maintainer='Jon Lund Steffensen',
    maintainer_email='jon_steffensen@uri.edu',
    url='https://github.com/zhanglab/psamm-import',
    license='GNU GPLv3+',

    classifiers = [
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
    ],

    packages=find_packages(),
    entry_points = {
        'console_scripts': [
            'psamm-import = psamm_import.importer:main',
            'psamm-import-bigg = psamm_import.bigg:main'
        ],
        'psamm.importer': [
            'iMA945 = psamm_import.datasource:ImportiMA945',
            'iRR1083 = psamm_import.datasource:ImportiRR1083',
            'iJO1366 = psamm_import.datasource:ImportiJO1366',
            'EColi_textbook = psamm_import.datasource:EColiTextbookImport',
            'STM_v1.0 = psamm_import.datasource:ImportSTMv1_0',
            'iJN746 = psamm_import.datasource:ImportiJN746',
            'iJP815 = psamm_import.datasource:ImportiJP815',
            'iSyn731 = psamm_import.datasource:ImportiSyn731',
            'iCce806 = psamm_import.datasource:ImportiCce806',
            'GSMN-TB = psamm_import.datasource:ImportGSMN_TB',
            'iNJ661 = psamm_import.datasource:ImportiNJ661',
            'iNJ661m = psamm_import.datasource:ImportiNJ661m',
            'iNJ661v = psamm_import.datasource:ImportiNJ661v',
            'iMR1_799 = psamm_import.datasource:ImportiMR1_799',
            'iMR4_812 = psamm_import.datasource:ImportiMR4_812',
            'iW3181_789 = psamm_import.datasource:ImportiW3181_789',
            'iOS217_672 = psamm_import.datasource:ImportiOS217_672',
            'ModelSEED = psamm_import.datasource:ImportModelSEED',
            'SBML-strict = psamm_import.datasource:SBMLStrictImporter',
            'SBML = psamm_import.datasource:SBMLNonstrictImporter',
            'COBRA-JSON = psamm_import.datasource:CobraJSONImporter'
        ]
    },

    install_requires=[
        'PyYAML>=3.11,<4.0',
        'xlrd',
        'psamm',
        'six'
    ])
