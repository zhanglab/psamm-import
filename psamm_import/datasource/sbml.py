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
# Copyright 2015-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""SBML importers."""

import os
import re
import glob
import logging

from six import iteritems, itervalues

from psamm.datasource import sbml
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry,
                                    DictReactionEntry as ReactionEntry)
from psamm.datasource.context import FilePathContext
from psamm.expression import boolean
from psamm.reaction import Compound, Reaction

from ..model import Importer, ParseError, ModelLoadError, MetabolicModel

logger = logging.getLogger(__name__)


class BaseImporter(Importer):
    """Base importer for reading metabolic model from an SBML file."""

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in SBML format.\n'
              'Expected files in source directory:\n'
              '- *.sbml')

    def _resolve_source(self, source):
        """Resolve source to filepath if it is a directory."""
        if os.path.isdir(source):
            sources = glob.glob(os.path.join(source, '*.sbml'))
            if len(sources) == 0:
                raise ModelLoadError('No .sbml file found in source directory')
            elif len(sources) > 1:
                raise ModelLoadError(
                    'More than one .sbml file found in source directory')
            return sources[0]
        return source

    def _get_reader(self, f):
        raise NotImplementedError('Subclasses must implement _get_reader()')

    def import_model(self, source):
        """Import and return model instance."""
        source = self._resolve_source(source)
        self._context = FilePathContext(source)
        with self._context.open() as f:
            self._reader = self._open_reader(f)

        model = MetabolicModel(self._reader.species, self._reader.reactions)
        model.compartments.update((c.id, c) for c in self._reader.compartments)

        if self._reader.name is not None:
            model.name = self._reader.name
        elif self._reader.id is not None:
            model.name = self._reader.id
        else:
            model.name = 'SBML'

        objective = self._reader.get_active_objective()
        if objective is not None:
            reactions = dict(objective.reactions)
            if len(reactions) != 1:
                logger.warning(
                    'Cannot convert objective {} since it does not'
                    ' consist of a single reaction'.format(objective.id))
            else:
                reaction, value = next(iteritems(reactions))
                if ((value > 0 and objective.type == 'minimize') or
                        (value < 0 and objective.type == 'maximize')):
                    logger.warning(
                        'Cannot convert objective {} since it is not a'
                        ' maximization objective'.format(objective.id))
                else:
                    logger.info('Detected biomass reaction: {}'.format(
                        reaction))
                    model.biomass_reaction = reaction

        return model


class StrictImporter(BaseImporter):
    """Read metabolic model from an SBML file using strict parser."""

    name = 'SBML-strict'
    title = 'SBML model (strict)'
    generic = True

    def _open_reader(self, f):
        try:
            return sbml.SBMLReader(f, strict=True, ignore_boundary=True)
        except sbml.ParseError as e:
            raise ParseError(e)


class NonstrictImporter(BaseImporter):
    """Read metabolic model from an SBML file using non-strict parser."""

    name = 'SBML'
    title = 'SBML model (non-strict)'
    generic = True

    _COBRA_ESCAPES = {
        '_DASH_': '-',
        '_FSLASH_': '/',
        '_BSLASH_': '\\',
        '_LPAREN_': '(',
        '_RPAREN_': ')',
        '_LSQBKT_': '[',
        '_RSQBKT_': ']',
        '_COMMA_': ',',
        '_PERIOD_': '.',
        '_APOS_': "'",
        '&amp;': '&',
        '&lt;': '<',
        '&gt;': '>',
        '&quot;': '"'
    }

    def _open_reader(self, f):
        try:
            return sbml.SBMLReader(f, strict=False, ignore_boundary=True)
        except sbml.ParseError as e:
            raise ParseError(e)

    def import_model(self, source):
        """Import and return model instance."""
        model = super(NonstrictImporter, self).import_model(source)

        objective_reactions = set()
        biomass_reaction = model.biomass_reaction
        for reaction in self._reader.reactions:
            # Check whether species occur multiple times
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

            if reaction.id not in model.limits:
                model.limits[reaction.id] = lower_bound, upper_bound
            else:
                current_lower, current_upper = model.limits[reaction.id]
                if current_lower is None:
                    current_lower = lower_bound
                if current_upper is None:
                    current_upper = upper_bound
                model.limits[reaction.id] = current_lower, current_upper

        if len(objective_reactions) == 1:
            biomass_reaction = next(iter(objective_reactions))
            logger.info('Detected biomass reaction: {}'.format(
                biomass_reaction))
        elif len(objective_reactions) > 1:
            logger.warning(
                'Multiple reactions are used as the'
                ' biomass reaction: {}'.format(objective_reactions))

        # Find compartments
        compartments = set()
        for reaction in itervalues(model.reactions):
            for compound, _ in reaction.equation.compounds:
                compartments.add(compound.compartment)

        # Detect prefixes
        compound_prefix = None
        if all(re.match(r'M_.+', c) for c in model.compounds):
            compound_prefix = 'M_'
            logger.info('Removing compound prefix {!r}'.format(
                compound_prefix))

        reaction_prefix = None
        if all(re.match(r'R_.+', r) for r in model.reactions):
            reaction_prefix = 'R_'
            logger.info('Removing reaction prefix {!r}'.format(
                reaction_prefix))

        compartment_prefix = None
        if all(re.match(r'C_.+', c) for c in compartments if c is not None):
            compartment_prefix = 'C_'
            logger.info('Removing compartment prefix {!r}'.format(
                compartment_prefix))

        compound_id_map = self._convert_compounds(
            model, compound_prefix=compound_prefix,
            compartment_prefix=compartment_prefix)
        reaction_id_map = self._convert_reactions(
            model, compound_id_map,
            reaction_prefix=reaction_prefix,
            compartment_prefix=compartment_prefix)

        # Remove prefix from biomass reaction
        if reaction_prefix is not None and model.biomass_reaction is not None:
            model.biomass_reaction = reaction_id_map.get(
                model.biomass_reaction, model.biomass_reaction)

        # Remove prefix from limits
        if reaction_prefix is not None:
            limits_to_check = set(model.limits)
            for reaction in limits_to_check:
                if reaction in reaction_id_map:
                    value = model.limits[reaction]
                    del model.limits[reaction]
                    new_reaction_id = reaction_id_map.get(reaction, reaction)
                    model.limits[new_reaction_id] = value

        return model

    def _parse_cobra_notes(self, element):
        for note in element.xml_notes.itertext():
            m = re.match(r'^([^:]+):(.+)$', note)
            if m:
                key, value = m.groups()
                key = key.strip().lower().replace(' ', '_')
                value = value.strip()
                m = re.match(r'^"(.*)"$', value)
                if m:
                    value = m.group(1)
                if value != '':
                    yield key, value

    def _convert_cobra_id(self, s):
        """Convert COBRA-specific symbol escapes in IDs."""
        for escape, symbol in iteritems(NonstrictImporter._COBRA_ESCAPES):
            s = s.replace(escape, symbol)
        return s

    def _convert_compounds(self, model, compound_prefix=None,
                           compartment_prefix=None):
        """Convert SBML species entries to compounds."""
        new_compounds = []
        id_map = {}
        for compound in itervalues(model.compounds):
            properties = dict(compound.properties)
            old_id = compound.id

            if compound_prefix is not None:
                if properties['id'].startswith(compound_prefix):
                    properties['id'] = properties['id'][len(compound_prefix):]

            if 'compartment' in properties and compartment_prefix is not None:
                if properties['compartment'].startswith(compartment_prefix):
                    properties['compartment'] = (
                        properties['compartment'][len(compartment_prefix):])

            properties['id'] = self._convert_cobra_id(properties['id'])

            # Extract information from notes
            if compound.xml_notes is not None:
                cobra_notes = dict(self._parse_cobra_notes(compound))

                for key in ('pubchem_id', 'chebi_id'):
                    if key in cobra_notes:
                        properties[key] = cobra_notes[key]

                if 'formula' in cobra_notes:
                    properties['formula'] = self._try_parse_formula(
                        compound.id, cobra_notes['formula'])

                if 'kegg_id' in cobra_notes:
                    properties['kegg'] = cobra_notes['kegg_id']

                if 'charge' in cobra_notes:
                    try:
                        value = int(cobra_notes['charge'])
                    except ValueError:
                        logger.warning(
                            'Unable to parse charge for {} as an'
                            ' integer: {}'.format(
                                compound.id, cobra_notes['charge']))
                    else:
                        properties['charge'] = value

            if old_id != properties['id']:
                id_map[old_id] = properties['id']
            new_compounds.append(CompoundEntry(
                properties, filemark=compound.filemark))

        model.compounds.clear()
        model.compounds.update((c.id, c) for c in new_compounds)

        return id_map

    def _convert_reactions(self, model, compound_id_map,
                           reaction_prefix=None, compartment_prefix=None):
        """Convert SBML reaction entries to reactions."""
        new_reactions = []
        id_map = {}
        for reaction in itervalues(model.reactions):
            properties = reaction.properties
            old_id = reaction.id

            if reaction_prefix is not None:
                if properties['id'].startswith(reaction_prefix):
                    properties['id'] = properties['id'][len(reaction_prefix):]

            properties['id'] = self._convert_cobra_id(properties['id'])

            # Extract information from notes
            if reaction.xml_notes is not None:
                cobra_notes = dict(self._parse_cobra_notes(reaction))

                if 'subsystem' in cobra_notes:
                    properties['subsystem'] = cobra_notes['subsystem']

                if 'gene_association' in cobra_notes:
                    assoc = self._try_parse_gene_association(
                        reaction.id, cobra_notes['gene_association'])
                    if assoc is not None:
                        properties['genes'] = assoc

                if 'ec_number' in cobra_notes:
                    properties['ec'] = cobra_notes['ec_number']

                if 'authors' in cobra_notes:
                    properties['authors'] = [
                        a.strip() for a in cobra_notes['authors'].split(';')]

                if 'confidence' in cobra_notes:
                    try:
                        value = int(cobra_notes['confidence'])
                    except ValueError:
                        logger.warning(
                            'Unable to parse confidence level for {} as an'
                            ' integer: {}'.format(
                                reaction.id, cobra_notes['confidence']))
                        value = cobra_notes['confidence']

                    properties['confidence'] = value

            # Convert compound IDs in reaction equation
            compounds = []
            for compound, value in properties['equation'].compounds:
                name = compound_id_map.get(compound.name, compound.name)

                compartment = compound.compartment
                if (compartment_prefix is not None and
                        compartment is not None and
                        compartment.startswith(compartment_prefix)):
                    compartment = compartment[len(compartment_prefix):]

                compounds.append(
                    (Compound(name, compartment=compartment), value))

                direction = properties['equation'].direction
                properties['equation'] = Reaction(direction, compounds)

            if old_id != properties['id']:
                id_map[old_id] = properties['id']
            new_reactions.append(ReactionEntry(
                properties, filemark=reaction.filemark))

        model.reactions.clear()
        model.reactions.update((r.id, r) for r in new_reactions)

        # Update model genes
        model.genes.clear()
        for reaction in itervalues(model.reactions):
            if reaction.genes is not None and isinstance(
                    reaction.genes, boolean.Expression):
                model.genes.update(v.symbol for v in reaction.genes.variables)

        return id_map
