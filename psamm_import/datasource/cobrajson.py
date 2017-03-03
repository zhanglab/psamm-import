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

"""Importer for the COBRApy JSON format."""

import os
import glob
import json
import logging
import decimal

from six import iteritems

from psamm.reaction import Reaction, Compound, Direction
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry,
                                    DictReactionEntry as ReactionEntry)

from ..model import Importer as BaseImporter, ModelLoadError, MetabolicModel

logger = logging.getLogger(__name__)


def _float_parser(num_str):
    num = float(num_str)
    if num.is_integer():
        return int(num)
    else:
        return decimal.Decimal(num_str)


class Importer(BaseImporter):
    """Read metabolic model from COBRApy JSON format."""

    name = 'json'
    title = 'COBRApy JSON'
    generic = True

    def help(self):
        """Print help text for importer."""
        print('Source must contain the model definition in COBRApy JSON'
              ' format.\n'
              'Expected files in source directory:\n'
              '- *.json')

    def _resolve_source(self, source):
        """Resolve source to filepath if it is a directory."""
        if os.path.isdir(source):
            sources = glob.glob(os.path.join(source, '*.json'))
            if len(sources) == 0:
                raise ModelLoadError('No .json file found in source directory')
            elif len(sources) > 1:
                raise ModelLoadError(
                    'More than one .json file found in source directory')
            return sources[0]
        return source

    def _import(self, file):
        model_doc = json.load(file, parse_float=_float_parser)
        model = MetabolicModel(
            self._read_compounds(model_doc), self._read_reactions(model_doc))
        model.name = model_doc.get('id', 'COBRA JSON model')

        biomass_reaction = None
        objective_reactions = set()
        for reaction in model_doc['reactions']:
            if reaction.get('objective_coefficient', 0) != 0:
                objective_reactions.add(reaction['id'])

        if len(objective_reactions) == 1:
            biomass_reaction = next(iter(objective_reactions))
            logger.info('Detected biomass reaction: {}'.format(
                biomass_reaction))
        elif len(objective_reactions) > 1:
            logger.warning(
                'Multiple reactions are used as the'
                ' biomass reaction: {}'.format(objective_reactions))

        model.biomass_reaction = biomass_reaction

        return model

    def _read_compounds(self, doc):
        for compound in doc['metabolites']:
            id = compound['id']
            name = compound.get('name')
            charge = compound.get('charge')

            formula_string = compound.get('formula')
            if formula_string is not None:
                formula = self._try_parse_formula(id, formula_string)
            else:
                formula = None

            yield CompoundEntry(dict(
                id=id, name=name, charge=charge, formula=formula))

    def _parse_reaction_equation(self, r):
        compounds = ((Compound(metabolite), value)
                     for metabolite, value in iteritems(r['metabolites']))
        if (r.get('lower_bound') == 0 and
                r.get('upper_bound') != 0):
            direction = Direction.Forward
        elif (r.get('lower_bound') != 0 and
              r.get('upper_bound') == 0):
            direction = Direction.Reverse
        else:
            direction = Direction.Both
        return Reaction(direction, compounds)

    def _read_reactions(self, doc):
        for reaction in doc['reactions']:
            id = reaction['id']
            name = reaction.get('name', None)
            equation = self._parse_reaction_equation(reaction)
            lower_flux = reaction.get('lower_bound')
            upper_flux = reaction.get('upper_bound')
            subsystem = reaction.get('subsystem')

            genes = reaction.get('gene_reaction_rule')
            if genes is not None:
                genes = self._try_parse_gene_association(id, genes)

            yield ReactionEntry(dict(
                id=id, name=name, equation=equation,
                lower_flux=lower_flux, upper_flux=upper_flux,
                subsystem=subsystem, genes=genes))

    def import_model(self, source):
        """Import and return model instance."""
        if not hasattr(source, 'read'):  # Not a File-like object
            with open(self._resolve_source(source), 'r') as f:
                return self._import(f)
        else:
            return self._import(source)
