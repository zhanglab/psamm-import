PSAMM model importers
=====================

PSAMM_ is an open source software that is designed for the curation and
analysis of metabolic models. It supports model version tracking, model
annotation, data integration, data parsing and formatting, consistency
checking, automatic gap filling, and model simulations.

The PSAMM model importers in this repository is collection of tools that allow
metabolic models to be converted from various formats to the internal
YAML-based format used by PSAMM.

The ``master`` branch tracks the latest release while the ``develop`` branch is
the latest version in development. Please apply any pull requests to the
``develop`` branch when creating the pull request.

.. _PSAMM: https://github.com/zhanglab/psamm

Overview
--------

This package provides two additional commands, ``psamm-import`` and
``psamm-import-bigg``. The ``psamm-import`` command allows the user to convert
various model files to YAML format. The conversion of SBML, COBRA JSON and
a number of Excel model formats are supported. The Excel format importers are
manually designed to load a specific model correctly, so only a limited set of
Excel models can be loaded. The SBML and COBRA JSON importers can be used with
any valid model in those formats. To see a list of all the supported importers,
use the following command:

.. code-block:: shell

    $ psamm-import list

To import an SBML model from the file ``ecoli_sbml_file.xml`` and extract it
as a YAML model in the directory ``ecoli_yaml``:

.. code-block:: shell

    $ psamm-import sbml --source ecoli_sbml_file.xml \
        --dest ecoli_yaml

The ``psamm-import-bigg`` command can be used to automatically download a model
from the BiGG_ online model database and convert it to YAML format. This
requires internet connection while running the command. To see a list of all
the available models in the database, use the following command:

.. code-block:: shell

    $ psamm-import-bigg list

To import the ``e_coli_core`` model and extract it as a YAML model in the
directory ``e_coli_core``:

.. code-block:: shell

    $ psamm-import-bigg e_coli_core --dest e_coli_core

.. _BiGG: http://bigg.ucsd.edu/

Install and documentation
-------------------------

Please see the main documentation for PSAMM that is available at
`Read the Docs`_.

.. _Read the Docs: https://psamm.readthedocs.org/

Software license
----------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

See LICENSE_.

.. _LICENSE: LICENSE
