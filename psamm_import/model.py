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

"""Common metabolic model representations for imported models.

This module contains various classes for representing the intermediate
result of parsing a model before it is converted to YAML format.
"""

from six import itervalues

from psamm.expression import boolean


class ImportError(Exception):
    """Exception used to signal a general import error."""


class ModelLoadError(ImportError):
    """Exception used to signal an error loading the model files."""


class ParseError(ImportError):
    """Exception used to signal an error parsing the model files."""


class _BaseEntry(object):
    """Base entry in loaded model."""

    def __init__(self, **kwargs):
        if 'id' not in kwargs:
            raise ValueError('No id was provided')
        self._values = kwargs
        self._id = self._values['id']

    @property
    def id(self):
        return self._id

    @property
    def properties(self):
        return dict(self._values)


class CompoundEntry(_BaseEntry):
    """Compound entry in loaded model."""

    @property
    def name(self):
        return self._values.get('name', None)

    @property
    def formula(self):
        return self._values.get('formula', None)

    @property
    def formula_neutral(self):
        return self._values.get('formula_neutral', None)

    @property
    def charge(self):
        return self._values.get('charge', None)

    @property
    def kegg(self):
        return self._values.get('kegg', None)

    @property
    def cas(self):
        return self._values.get('cas', None)


class ReactionEntry(_BaseEntry):
    """Reaction entry in loaded model."""

    @property
    def name(self):
        return self._values.get('name', None)

    @property
    def genes(self):
        return self._values.get('genes', None)

    @property
    def equation(self):
        return self._values.get('equation', None)

    @property
    def subsystem(self):
        return self._values.get('subsystem', None)

    @property
    def ec(self):
        return self._values.get('ec', None)


class MetabolicModel(object):
    def __init__(self, name, compounds, reactions):
        self._name = name
        self._compounds = dict((c.id, c) for c in compounds)
        self._reactions = dict((r.id, r) for r in reactions)
        self._biomass_reaction = None

        self._genes = set()
        for r in itervalues(self._reactions):
            if hasattr(r, 'genes') and r.genes is not None:
                if isinstance(r.genes, boolean.Expression):
                    self._genes.update(g.symbol for g in r.genes.variables)
                else:
                    self._genes.update(r.genes)

        self._check_reaction_compounds()

    @property
    def name(self):
        return self._name

    @property
    def reactions(self):
        return self._reactions

    @property
    def compounds(self):
        return self._compounds

    @property
    def genes(self):
        return self._genes

    @property
    def biomass_reaction(self):
        return self._biomass_reaction

    @biomass_reaction.setter
    def biomass_reaction(self, value):
        if value is not None and value not in self._reactions:
            raise ValueError('Invalid reaction')
        self._biomass_reaction = value

    def print_summary(self):
        """Print model summary."""
        print('Model: {}'.format(self.name))
        print('- Biomass reaction: {}'.format(self.biomass_reaction))
        print('- Compounds: {}'.format(len(self.compounds)))
        print('- Reactions: {}'.format(len(self.reactions)))
        print('- Genes: {}'.format(len(self.genes)))

    def _check_reaction_compounds(self):
        """Check that reaction compounds are defined in the model."""
        undefined = set()
        for reaction in itervalues(self.reactions):
            if reaction.equation is not None:
                for compound, value in reaction.equation.compounds:
                    if compound.name not in self.compounds:
                        undefined.add((reaction.id, compound.name))

        if len(undefined) > 0:
            raise ParseError('Some reaction compounds are not defined in'
                             ' the model: {}'.format(undefined))


class Importer(object):
    """Base importer class."""
