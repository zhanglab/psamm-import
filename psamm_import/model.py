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
# Copyright 2015-2016  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Common metabolic model representations for imported models.

This module contains various classes for representing the intermediate
result of parsing a model before it is converted to YAML format.
"""

import logging
from collections import Counter

from six import iteritems, itervalues, text_type

from psamm.expression import boolean
from psamm.datasource.reaction import (parse_reaction,
                                       ParseError as ReactionParseError)
from psamm import formula
from psamm.reaction import Direction

logger = logging.getLogger(__name__)


class ImportError(Exception):
    """Exception used to signal a general import error."""


class ModelLoadError(ImportError):
    """Exception used to signal an error loading the model files."""


class ParseError(ImportError):
    """Exception used to signal an error parsing the model files."""


class _BaseEntry(object):
    """Base entry in loaded model."""

    def __init__(self, **kwargs):
        self._values = {key: value for key, value in iteritems(kwargs)
                        if value is not None}
        if 'id' not in self._values:
            raise ValueError('No id was provided')
        self._id = self._values['id']

    @property
    def id(self):
        return self._id

    @property
    def properties(self):
        return self._values


class CompoundEntry(_BaseEntry):
    """Compound entry in loaded model."""

    @property
    def name(self):
        """Compound name."""
        return self._values.get('name', None)

    @property
    def formula(self):
        """Compound formula."""
        return self._values.get('formula', None)

    @property
    def formula_neutral(self):
        """Compound formula (neutral)."""
        return self._values.get('formula_neutral', None)

    @property
    def charge(self):
        """Compound charge."""
        return self._values.get('charge', None)

    @property
    def kegg(self):
        """KEGG identifier."""
        return self._values.get('kegg', None)

    @property
    def cas(self):
        """CAS identifier."""
        return self._values.get('cas', None)


class ReactionEntry(_BaseEntry):
    """Reaction entry in loaded model."""

    @property
    def name(self):
        """Reaction name."""
        return self._values.get('name', None)

    @property
    def genes(self):
        """Gene list or boolean expression."""
        return self._values.get('genes', None)

    @property
    def equation(self):
        """Reaction equation."""
        return self._values.get('equation', None)

    @property
    def subsystem(self):
        """Reaction subsystem classification."""
        return self._values.get('subsystem', None)

    @property
    def ec(self):
        """EC classifier."""
        return self._values.get('ec', None)


class MetabolicModel(object):
    """Intermediate model representation parsed from external source."""

    def __init__(self, compounds, reactions):
        """Create model with name, compounds and reactions."""
        self._compounds = dict((c.id, c) for c in compounds)

        self._reactions = {}
        self._limits = {}
        for reaction in reactions:
            # Extract flux limits embedded in reaction properties
            lower_flux = reaction.properties.pop('lower_flux', None)
            upper_flux = reaction.properties.pop('upper_flux', None)
            if lower_flux is not None or upper_flux is not None:
                self._limits[reaction.id] = lower_flux, upper_flux

            self._reactions[reaction.id] = reaction

        self._medium = {}

        self._name = None
        self._biomass_reaction = None
        self._extracellular_compartment = None
        self._default_compartment = None
        self._default_flux_limit = None

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
        """Model name."""
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def reactions(self):
        """Model reactions as dictionary."""
        return self._reactions

    @property
    def compounds(self):
        """Model compounds as dictionary."""
        return self._compounds

    @property
    def genes(self):
        """Model genes as set."""
        return self._genes

    @property
    def limits(self):
        """Return limits on reaction fluxes."""
        return self._limits

    @property
    def medium(self):
        """Medium definition."""
        return self._medium

    @property
    def biomass_reaction(self):
        """Main biomass reaction of model."""
        return self._biomass_reaction

    @biomass_reaction.setter
    def biomass_reaction(self, value):
        if value is not None and value not in self._reactions:
            raise ValueError('Invalid reaction')
        self._biomass_reaction = value

    @property
    def extracellular_compartment(self):
        """Extracellular compartment specified by the model."""
        return self._extracellular_compartment

    @extracellular_compartment.setter
    def extracellular_compartment(self, value):
        self._extracellular_compartment = value

    @property
    def default_compartment(self):
        """Default compartment specified by the model."""
        return self._default_compartment

    @default_compartment.setter
    def default_compartment(self, value):
        self._default_compartment = value

    @property
    def default_flux_limit(self):
        """Default flux limit specified by the model."""
        return self._default_flux_limit

    @default_flux_limit.setter
    def default_flux_limit(self, value):
        self._default_flux_limit = value

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
        for reaction in itervalues(self._reactions):
            if reaction.equation is not None:
                for compound, value in reaction.equation.compounds:
                    if compound.name not in self._compounds:
                        undefined.add((reaction.id, compound.name))

        if len(undefined) > 0:
            raise ParseError('Some reaction compounds are not defined in'
                             ' the model: {}'.format(undefined))


class Importer(object):
    """Base importer class."""

    def _try_parse_formula(self, compound_id, s):
        """Try to parse the given compound formula string.

        Logs a warning if the formula could not be parsed.
        """
        s = s.strip()
        if s == '':
            return None

        try:
            # Do not return the parsed formula. For now it is better to keep
            # the original formula string unchanged in all cases.
            formula.Formula.parse(s)
        except ValueError:
            logger.warning('Unable to parse compound formula {}: {}'.format(
                compound_id, s))

        return s

    def _try_parse_reaction(self, reaction_id, s,
                            parser=parse_reaction, **kwargs):
        """Try to parse the given reaction equation string.

        Returns the parsed Reaction object, or raises an error if the reaction
        could not be parsed.
        """
        try:
            return parser(s, **kwargs)
        except ReactionParseError as e:
            if e.indicator is not None:
                logger.error(u'{}\n{}\n{}'.format(
                    str(e), s, e.indicator))
            raise ParseError('Unable to parse reaction {}: {}'.format(
                reaction_id, s))

    def _try_parse_gene_association(self, reaction_id, s):
        """Try to parse the given gene association rule.

        Logs a warning if the association rule could not be parsed and returns
        the original string. Otherwise, returns the boolean.Expression object.
        """
        if s == '':
            return None

        try:
            return boolean.Expression(s)
        except boolean.ParseError as e:
            msg = u'Failed to parse gene association for {}: {}'.format(
                reaction_id, text_type(e))
            if e.indicator is not None:
                msg += u'\n{}\n{}'.format(s, e.indicator)
            logger.warning(msg)

        return s


def detect_extracellular_compartment(model):
    """Detect the identifier for equations with extracellular compartments."""
    extracellular_key = Counter()

    for reaction_id, reaction in iteritems(model.reactions):
        if 'equation' not in reaction.properties:
            continue

        equation = reaction.properties['equation']
        if len(equation.compounds) == 1:
            compound, _ = equation.compounds[0]
            compartment = compound.compartment
            extracellular_key[compartment] += 1
    if len(extracellular_key) == 0:
        return None
    else:
        best_key, _ = extracellular_key.most_common(1)[0]
    logger.info('{} is extracellular compartment'.format(best_key))
    return best_key


def convert_exchange_to_medium(model):
    """Convert exchange reactions in model to medium definition.

    Only exchange reactions in the extracellular compartment are converted to
    medium. The extracelluar compartment must be defined for the model.
    """
    # Build set of exchange reactions
    exchanges = set()
    for reaction_id, reaction in iteritems(model.reactions):
        equation = reaction.properties.get('equation')
        if equation is None:
            continue

        if len(equation.compounds) != 1:
            # Provide warning for exchange reactions with more than
            # one compound, they won't be put into the medium definition
            if (len(equation.left) == 0) != (len(equation.right) == 0):
                logger.warning('Exchange reaction {} has more than one'
                               ' compound, it was not converted to'
                               ' medium compounds'.format(reaction.id))
            continue

        exchanges.add(reaction_id)

    # Convert exchange reactions into medium
    for reaction_id in exchanges:
        equation = model.reactions[reaction_id].properties['equation']
        compound, value = equation.compounds[0]
        if compound.compartment != model.extracellular_compartment:
            continue

        if compound in model.medium:
            logger.warning(
                'Compound {} is already defined in the medium'.format(
                    compound))
            continue

        # We multiply the flux bounds by value in order to create equivalent
        # exchange reactions with stoichiometric value of one. If the flux
        # bounds are not set but the reaction is unidirectional, the implicit
        # flux bounds must be used.
        lower_flux, upper_flux = None, None
        if reaction_id in model.limits:
            lower, upper = model.limits[reaction_id]
            if lower is not None:
                lower_flux = lower * abs(value)
            if upper is not None:
                upper_flux = upper * abs(value)
        elif equation.direction == Direction.Forward:
            lower_flux = 0.0
        elif equation.direction == Direction.Reverse:
            upper_flux = 0.0

        # If the stoichiometric value of the reaction is reversed, the flux
        # limits must be flipped.
        if value > 0:
            lower_flux, upper_flux = (
                -upper_flux if upper_flux is not None else None,
                -lower_flux if lower_flux is not None else None)

        model.medium[compound] = reaction_id, lower_flux, upper_flux

        del model.reactions[reaction_id]
        model.limits.pop(reaction_id, None)
