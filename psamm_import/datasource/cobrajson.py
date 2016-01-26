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

"""Importer for the COBRApy JSON format."""

import os
import glob
import json
import logging

from six import iteritems

from psamm.reaction import Reaction, Compound, Direction

from ..model import (Importer as BaseImporter, ModelLoadError,
                     ParseError, CompoundEntry, ReactionEntry, MetabolicModel)

logger = logging.getLogger(__name__)


class Importer(BaseImporter):
    """Read metabolic model from COBRApy JSON format"""

    name = 'json'
    title = 'COBRApy JSON'
    generic = True

    def help(self):
        print('Source must contain the model definition in COBRApy JSON'
              ' format.\n'
              'Expected files in source directory:\n'
              '- *.json')

    def _resolve_source(self, source):
        """Resolve source to filepath if it is a directory"""
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
        model_doc = json.load(file)
        model = MetabolicModel(
            model_doc.get('id', 'COBRA JSON model'),
            self._read_compounds(model_doc), self._read_reactions(model_doc))

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
            formula = self._try_parse_formula(id, compound.get('formula'))
            yield CompoundEntry(id=id, name=name, charge=charge,
                                formula=formula)

    def _parse_reaction_equation(self, doc):
        compounds = ((Compound(metabolite), value)
                     for metabolite, value in iteritems(doc))
        return Reaction(Direction.Both, compounds)


    def _read_reactions(self, doc):
        for reaction in doc['reactions']:
            id = reaction['id']
            name = reaction.get('name', None)
            equation = self._parse_reaction_equation(reaction['metabolites'])
            lower_flux = reaction.get('lower_bound')
            upper_flux = reaction.get('upper_bound')
            subsystem = reaction.get('subsystem')

            genes = reaction.get('gene_reaction_rule')
            if genes is not None:
                genes = self._try_parse_gene_association(id, genes)

            yield ReactionEntry(id=id, name=name, equation=equation,
                                lower_flux=lower_flux, upper_flux=upper_flux,
                                subsystem=subsystem, genes=genes)

    def import_model(self, source):
        if not hasattr(source, 'read'):  # Not a File-like object
            with open(self._resolve_source(source), 'r') as f:
                return self._import(f)
        else:
            return self._import(source)
