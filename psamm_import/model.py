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

"""Common metabolic model representations for imported models.

This module contains various classes for representing the intermediate
result of parsing a model before it is converted to YAML format.
"""

import logging
from collections import Counter
from itertools import product

from six import iteritems, itervalues, text_type

from psamm.expression import boolean
from psamm.datasource.reaction import (parse_reaction,
                                       ParseError as ReactionParseError)
from psamm.datasource.entry import DictCompoundEntry, DictCompartmentEntry
from psamm import formula
from psamm.reaction import Direction

logger = logging.getLogger(__name__)


class ImportError(Exception):
    """Exception used to signal a general import error."""


class ModelLoadError(ImportError):
    """Exception used to signal an error loading the model files."""


class ParseError(ImportError):
    """Exception used to signal an error parsing the model files."""


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

        self._exchange = {}

        self._compartments = {}
        self._compartment_adjacency = {}

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
    def compartments(self):
        """Model compartments as dictionary."""
        return self._compartments

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
    def exchange(self):
        """Return exchange definition."""
        return self._exchange

    @property
    def compartment_adjacency(self):
        """Return compartment adjacency dictionary."""
        return self._compartment_adjacency

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
        print('- Compartments: {}'.format(len(self.compartments)))
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
        equation = reaction.equation
        if equation is None:
            continue

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


def convert_exchange_to_compounds(model):
    """Convert exchange reactions in model to exchange compounds.

    Only exchange reactions in the extracellular compartment are converted.
    The extracelluar compartment must be defined for the model.
    """
    # Build set of exchange reactions
    exchanges = set()
    for reaction_id, reaction in iteritems(model.reactions):
        equation = reaction.properties.get('equation')
        if equation is None:
            continue

        if len(equation.compounds) != 1:
            # Provide warning for exchange reactions with more than
            # one compound, they won't be put into the exchange definition
            if (len(equation.left) == 0) != (len(equation.right) == 0):
                logger.warning('Exchange reaction {} has more than one'
                               ' compound, it was not converted to'
                               ' exchange compound'.format(reaction.id))
            continue

        exchanges.add(reaction_id)

    # Convert exchange reactions into exchange compounds
    for reaction_id in exchanges:
        equation = model.reactions[reaction_id].properties['equation']
        compound, value = equation.compounds[0]
        if compound.compartment != model.extracellular_compartment:
            continue

        if compound in model.exchange:
            logger.warning(
                'Compound {} is already defined in the exchange'
                ' definition'.format(compound))
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

        if lower_flux is None and equation.direction == Direction.Forward:
            lower_flux = 0.0
        if upper_flux is None and equation.direction == Direction.Reverse:
            upper_flux = 0.0

        # If the stoichiometric value of the reaction is reversed, the flux
        # limits must be flipped.
        if value > 0:
            lower_flux, upper_flux = (
                -upper_flux if upper_flux is not None else None,
                -lower_flux if lower_flux is not None else None)

        model.exchange[compound] = reaction_id, lower_flux, upper_flux

        del model.reactions[reaction_id]
        model.limits.pop(reaction_id, None)


def infer_compartment_entries(model):
    """Infer compartment entries for model based on reaction compounds."""
    compartment_ids = set()
    for reaction in itervalues(model.reactions):
        equation = reaction.equation
        if equation is None:
            continue

        for compound, _ in equation.compounds:
            compartment = compound.compartment
            if compartment is None:
                compartment = model.default_compartment

            if compartment is not None:
                compartment_ids.add(compartment)

    for compartment in compartment_ids:
        if compartment in model.compartments:
            continue

        entry = DictCompartmentEntry(dict(id=compartment))
        model.compartments[entry.id] = entry


def infer_compartment_adjacency(model):
    """Infer compartment adjacency for model based on reactions."""
    def reaction_compartments(seq):
        for compound, _ in seq:
            compartment = compound.compartment
            if compartment is None:
                compartment = model.default_compartment

            if compartment is not None:
                yield compartment

    for reaction in itervalues(model.reactions):
        equation = reaction.equation
        if equation is None:
            continue

        left = reaction_compartments(equation.left)
        right = reaction_compartments(equation.right)
        for c1, c2 in product(left, right):
            if c1 == c2:
                continue
            model.compartment_adjacency.setdefault(c1, set()).add(c2)
            model.compartment_adjacency.setdefault(c2, set()).add(c1)


def merge_equivalent_compounds(model):
    """Merge equivalent compounds in various compartments.

    Tries to detect and merge compound entries that represent the same
    compound in different compartments. The entries are only merged if all
    properties are equivalent. Compound entries must have an ID with a suffix
    of an underscore followed by the compartment ID. This suffix will be
    stripped and compounds with identical IDs are merged if the properties
    are identical.
    """
    def dicts_are_compatible(d1, d2):
        return all(key not in d1 or key not in d2 or d1[key] == d2[key]
                   for key in set(d1) | set(d2))

    compound_compartment = {}
    inelegible = set()
    for reaction in itervalues(model.reactions):
        equation = reaction.equation
        if equation is None:
            continue

        for compound, _ in equation.compounds:
            compartment = compound.compartment
            if compartment is not None:
                compound_compartment[compound.name] = compartment
                if not compound.name.endswith('_{}'.format(compartment)):
                    inelegible.add(compound.name)

    compound_groups = {}
    for compound_id, compartment in iteritems(compound_compartment):
        if compound_id in inelegible:
            continue

        suffix = '_{}'.format(compound_compartment[compound_id])
        if compound_id.endswith(suffix):
            group_name = compound_id[:-len(suffix)]
            compound_groups.setdefault(group_name, set()).add(compound_id)

    compound_mapping = {}
    merged_compounds = {}
    for group, compound_set in iteritems(compound_groups):
        # Try to merge as many compounds as possible
        merged = []
        for compound_id in compound_set:
            props = dict(model.compounds[compound_id].properties)
            props.pop('id', None)
            props.pop('compartment', None)
            for merged_props, merged_set in merged:
                if dicts_are_compatible(props, merged_props):
                    merged_set.add(compound_id)
                    merged_props.update(props)
                    break
                else:
                    keys = set(key for key in set(props) | set(merged_props)
                               if key not in props or
                               key not in merged_props or
                               props[key] != merged_props[key])
                    logger.info(
                        'Unable to merge {} into {}, difference in'
                        ' keys: {}'.format(
                            compound_id, ', '.join(merged_set),
                            ', '.join(keys)))
            else:
                merged.append((props, {compound_id}))

        if len(merged) == 1:
            # Merge into one set with the group name
            merged_props, merged_set = merged[0]

            for compound_id in merged_set:
                compound_mapping[compound_id] = group
            merged_compounds[group] = merged_props
        else:
            # Since we cannot merge all compounds, create new group names
            # based on the group and compartments.
            for merged_props, merged_set in merged:
                compartments = set(compound_compartment[c] for c in merged_set)
                merged_name = '{}_{}'.format(
                    group, '_'.join(sorted(compartments)))

                for compound_id in merged_set:
                    compound_mapping[compound_id] = merged_name
                merged_compounds[merged_name] = merged_props

    # Translate reaction compounds
    for reaction in itervalues(model.reactions):
        equation = reaction.equation
        if equation is None:
            continue

        reaction.properties['equation'] = equation.translated_compounds(
            lambda c: compound_mapping.get(c, c))

    # Translate compound entries
    new_compounds = []
    for compound in itervalues(model.compounds):
        if compound.id not in compound_mapping:
            new_compounds.append(compound)
        else:
            group = compound_mapping[compound.id]
            if group not in merged_compounds:
                continue
            props = merged_compounds.pop(group)
            props['id'] = group
            new_compounds.append(DictCompoundEntry(
                props, filemark=compound.filemark))

    model.compounds.clear()
    model.compounds.update((c.id, c) for c in new_compounds)

    # Translate exchange
    for compound in compound_mapping:
        if compound in model.exchange:
            value = model.exchange.pop(compound)
            model.exchange[compound_mapping[compound]] = value
