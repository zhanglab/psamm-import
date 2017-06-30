PSAMM Excel model importers
===========================

**Note: This package no longer includes the SBML and COBRA JSON importers.**
The SBML and COBRA JSON importers have been merged into the main PSAMM
package and the import commands (``psamm-import`` and ``psamm-import-bigg``)
are now available from the PSAMM package.

PSAMM_ is an open source software that is designed for the curation and
analysis of metabolic models. It supports model version tracking, model
annotation, data integration, data parsing and formatting, consistency
checking, automatic gap filling, and model simulations.

The PSAMM Excel model importers in this repository is a collection of tools
that allow metabolic models to be converted from model-specific Excel formats
to the internal YAML-based format used by PSAMM.

The ``master`` branch tracks the latest release while the ``develop`` branch is
the latest version in development. Please apply any pull requests to the
``develop`` branch when creating the pull request.

.. _PSAMM: https://github.com/zhanglab/psamm

Overview
--------

This package provides additional Excel format importers that are manually
designed to load a specific model correctly, so only a limited set of Excel
models can be loaded. After installing this package, the Excel importers will
become available from ``psamm-import``. Use the following command to see a
list of models that can be imported:

.. code-block:: shell

    $ psamm-import list

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
