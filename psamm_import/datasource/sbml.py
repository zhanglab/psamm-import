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

"""SBML importers."""

import os
import re
import glob
import logging

from psamm.datasource import sbml

from ..model import (Importer, ParseError, CompoundEntry, ReactionEntry,
                     MetabolicModel)

logger = logging.getLogger(__name__)


class BaseImporter(Importer):
    """Base importer for reading metabolic model from an SBML file."""

    def help(self):
        print('Source must contain the model definition in SBML format.\n'
              'Expected files in source directory:\n'
              '- *.sbml')

    def _resolve_source(self, source):
        """Resolve source to filepath if it is a directory"""
        if os.path.isdir(source):
            sources = glob.glob(os.path.join(source, '*.sbml'))
            if len(sources) == 0:
                raise ParseError('No .sbml file found in source directory')
            elif len(sources) > 1:
                raise ParseError(
                    'More than one .sbml file found in source directory')
            return sources[0]
        return source

    def _get_reader(self, f):
        raise NotImplementedError('Subclasses must implement _get_reader()')

    def import_model(self, source):
        source = self._resolve_source(source)
        with open(source, 'r') as f:
            self._reader = self._open_reader(f)

        model_name = 'SBML'
        if self._reader.name is not None:
            model_name = self._reader.name
        elif self._reader.id is not None:
            model_name = self._reader.id

        model = MetabolicModel(
            model_name, self._reader.species, self._reader.reactions)

        return model


class StrictImporter(BaseImporter):
    """Read metabolic model from an SBML file using strict parser."""

    name = 'SBML-strict'
    title = 'SBML model (strict)'

    def _open_reader(self, f):
        return sbml.SBMLReader(f, strict=True, ignore_boundary=True)


class NonstrictImporter(BaseImporter):
    """Read metabolic model from an SBML file using non-strict parser."""

    name = 'SBML'
    title = 'SBML model (non-strict)'

    def _open_reader(self, f):
        return sbml.SBMLReader(f, strict=False, ignore_boundary=True)

    def import_model(self, source):
        model = super(NonstrictImporter, self).import_model(source)

        biomass_reaction = None
        objective_reactions = set()
        flux_limits = {}
        for reaction in self._reader.reactions:
            # Check whether species multiple times
            compounds = set()
            for c, _ in reaction.equation.compounds:
                if c.name in compounds:
                    logger.warning(
                        'Compound {} appears multiple times in the same'
                        ' reaction: {}'.format(c.name, reaction.id))
                compounds.add(c.name)

            # Try to detect biomass reaction and ambiguous flux bounds
            upper_bound = None
            lower_bound = None
            for parameter in reaction.kinetic_law_reaction_parameters:
                pid, name, value, units = parameter
                if (pid == 'OBJECTIVE_COEFFICIENT' or
                        name == 'OBJECTIVE_COEFFICIENT'):
                    try:
                        value = float(value)
                    except ValueError:
                        continue
                    if value != 0:
                        objective_reactions.add(reaction.id)
                elif pid == 'UPPER_BOUND' or name == 'UPPER_BOUND':
                    try:
                        value = float(value)
                    except ValueError:
                        continue
                    upper_bound = value
                elif pid == 'LOWER_BOUND' or name == 'LOWER_BOUND':
                    try:
                        value = float(value)
                    except ValueError:
                        continue
                    lower_bound = value

            if upper_bound is not None and upper_bound <= 0:
                logger.warning('Upper bound of {} is {}'.format(
                    reaction.id, upper_bound))
            if (lower_bound is not None and
                    lower_bound < 0 and not reaction.reversible):
                logger.warning('Lower bound of irreversible reaction {} is'
                               ' {}'.format(reaction.id, lower_bound))

            flux_limits[reaction.id] = (lower_bound, upper_bound)

        if len(objective_reactions) == 1:
            biomass_reaction = next(iter(objective_reactions))
            logger.info('Detected biomass reaction: {}'.format(
                biomass_reaction))
        elif len(objective_reactions) > 1:
            logger.warning(
                'Multiple reactions are used as the'
                ' biomass reaction: {}'.format(objective_reactions))

        model = MetabolicModel(
            model.name,
            self._convert_compounds(model.compounds.itervalues()),
            self._convert_reactions(model.reactions.itervalues(),
                                    flux_limits))
        model.biomass_reaction = biomass_reaction

        return model

    def _convert_compounds(self, compounds):
        """Convert SBML species entries to compounds."""
        for compound in compounds:
            properties = compound.properties

            # Extract information from notes
            if compound.xml_notes is not None:
                for note in compound.xml_notes.itertext():
                    m = re.match(r'FORMULA: (.+)$', note)
                    if m:
                        properties['formula'] = m.group(1)

                    m = re.match(r'KEGG ID: (.+)$', note)
                    if m:
                        properties['kegg'] = m.group(1)

                    m = re.match(r'PubChem ID: (.+)$', note)
                    if m:
                        properties['pubchem_id'] = m.group(1)

                    m = re.match(r'ChEBI ID: (.+)$', note)
                    if m:
                        properties['chebi_id'] = m.group(1)

            yield CompoundEntry(**properties)

    def _convert_reactions(self, reactions, flux_limits):
        """Convert SBML reaction entries to reactions."""
        for reaction in reactions:
            properties = reaction.properties

            # Extract information from notes
            if reaction.xml_notes is not None:
                for note in reaction.xml_notes.itertext():
                    m = re.match(r'SUBSYSTEM: (.+)$', note)
                    if m:
                        properties['subsystem'] = m.group(1)

                    m = re.match(r'GENE_ASSOCIATION: (.+)$', note)
                    if m:
                        properties['gene_association'] = m.group(1)

            # Extract flux limits
            if reaction.id in flux_limits:
                lower, upper = flux_limits[reaction.id]
                if lower is not None:
                    properties['lower_flux'] = lower
                if upper is not None:
                    properties['upper_flux'] = upper

            yield ReactionEntry(**properties)
