
"""Functions related to loading models"""

from __future__ import print_function

import os
import glob
import re
import csv
from collections import namedtuple
import logging

import xlrd
from psamm.datasource.misc import (parse_metnet_reaction,
                                   parse_sudensimple_reaction)
from psamm.datasource import modelseed, sbml
from psamm.formula import Formula
from psamm.reaction import Reaction, Compound

logger = logging.getLogger(__name__)


class _BaseEntry(object):
    """Base entry in loaded model"""

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
    """Compound entry in loaded model"""

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
    """Reaction entry in loaded model"""

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


class Importer(object):
    """Base importer class"""


class MetabolicModel(object):
    def __init__(self, name, compounds, reactions):
        self._name = name
        self._compounds = dict((c.id, c) for c in compounds)
        self._reactions = dict((r.id, r) for r in reactions)
        self._biomass_reaction = None

        self._genes = set()
        for r in self._reactions.itervalues():
            if hasattr(r, 'genes') and r.genes is not None:
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
        """Print model summary"""
        print('Model: {}'.format(self.name))
        print('- Biomass reaction: {}'.format(self.biomass_reaction))
        print('- Compounds: {}'.format(len(self.compounds)))
        print('- Reactions: {}'.format(len(self.reactions)))
        print('- Genes: {}'.format(len(self.genes)))

    def _check_reaction_compounds(self):
        """Check that reaction compounds are defined in the model"""
        undefined = set()
        for reaction in self.reactions.itervalues():
            if reaction.equation is not None:
                for compound, value in reaction.equation.compounds:
                    if compound.name not in self.compounds:
                        undefined.add((reaction.id, compound.name))

        if len(undefined) > 0:
            raise ParseError('Some reaction compounds are not defined in'
                             'the model: {}'.format(undefined))


class ParseError(Exception):
    """Exception used to signal a parsing error"""


class ImportiMA945(Importer):
    name = 'iMA945'
    title = 'Salmonella enterica iMA945 (Excel format), AbuOun et al., 2009'

    filename = 'jbc.M109.005868-5.xls'

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())
        model.biomass_reaction = 'ST_biomass_core'

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('compounds')
        for i in range(1, sheet.nrows):
            compound_id, name, formula, charge, cas, formula_neutral, kegg = (
                sheet.row_values(i))

            if compound_id.strip() == '':
                continue

            # Fixup model errors
            m = re.match(r'^(.*)(C\d{5})$', formula_neutral)
            if m:
                formula_neutral = m.group(1)
                kegg = m.group(2)

            name = None if name == '' else name
            name = re.sub(r'\'', '', name)
            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            if formula_neutral.strip() != '':
                formula_neutral = Formula.parse(formula_neutral)
                formula = formula_neutral
            else:
                formula_neutral = None
                if formula.strip() != '':
                    m = re.match(r'^(.*)-\d$', formula)
                    if m:
                        formula = m.group(1)
                    formula = Formula.parse(formula)
                else:
                    formula = None

            kegg = None if kegg == '' else kegg

            yield CompoundEntry(id=compound_id, name=name,
                                formula=formula,
                                formula_neutral=formula_neutral,
                                charge=charge, kegg=kegg)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('reactions')
        for i in range(1, sheet.nrows):
            reaction_id, name, equation, genes = (
                sheet.row_values(i, end_colx=4))

            if reaction_id.strip() == '':
                continue

            # Fixup model errors
            if reaction_id in ('FACOAL100t2pp', 'FACOAL80t2pp'):
                # Reaction equation and genes are messed up
                equation = re.sub(r'STM1818$', ']', equation)
                genes = 'STM1818'
            elif reaction_id == 'FACOAL60t2pp':
                equation = re.sub(r'STM1818$', '', equation)
                genes = 'STM1818'
            elif reaction_id == 'NTRIR4pp':
                m = re.match(r'^(.*nh4\[p\])(.*)$', equation)
                equation = m.group(1)
                genes = m.group(2)
            elif reaction_id in ('FE3DHBZSabcpp', '14GLUCANabcpp'):
                m = re.match(r'^(.*pi\[c\])(.*)$', equation)
                equation = m.group(1)
                genes = m.group(2)

            name = None if name == '' else name

            if equation != '':
                equation = parse_metnet_reaction(equation)
            else:
                equation = None

            if genes != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'STM\d+', genes))
            else:
                genes = None

            yield ReactionEntry(id=reaction_id, name=name,
                                genes=genes, equation=equation)


class ImportiRR1083(Importer):
    name = 'iRR1083'
    title = ('Salmonella enterica iRR1083 (Excel format),'
             ' Raghunathan et al., 2009')

    filename = '1752-0509-3-38-s1.xls'

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('Metabolites')
        for i in range(1, sheet.nrows):
            compound_id, name, formula_neutral, charge, kegg = (
                sheet.row_values(i))

            if compound_id.strip() == '':
                continue

            name = None if name == '' else name

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            if formula_neutral.strip() != '':
                formula_neutral = Formula.parse(formula_neutral)
            else:
                formula_neutral = None

            kegg = None if kegg == '' else kegg

            yield CompoundEntry(id=compound_id, name=name,
                                formula=formula_neutral,
                                formula_neutral=formula_neutral,
                                charge=charge, kegg=kegg)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('Gene Protein Reaction iRR1083')
        for i in range(3, sheet.nrows):
            genes, protein, reaction_id, name, equation, subsystem = (
                sheet.row_values(i, end_colx=6))

            if reaction_id.strip() == '':
                continue

            if genes != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'STM\d+', genes))
            else:
                genes = None

            protein = None if protein == '' else protein
            name = None if name == '' else name

            if equation != '':
                equation = parse_metnet_reaction(equation)
            else:
                equation = None

            subsystem = None if subsystem == '' else subsystem

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation)


class ImportiJO1366(Importer):
    name = 'iJO1366'
    title = ('Escerichia coli iJO1366 (Excel format),'
             ' Orth et al., 2011')

    filename = 'inline-supplementary-material-2.xls'

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())
        model.biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M'

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('Table 3')
        for i in range(1, sheet.nrows):
            (compound_id, name, formula_neutral, formula, charge, compartment,
                kegg, cas, alt_names) = sheet.row_values(i, end_colx=9)

            if compound_id.strip() == '':
                continue

            compound_id = re.match(r'^(.*)\[.\]$', compound_id).group(1)
            name = None if name.strip() == '' else name

            if formula_neutral.strip() != '':
                formula_neutral = Formula.parse(formula_neutral)
                formula = formula_neutral
            else:
                formula_neutral = None
                if formula.strip() != '':
                    formula = Formula.parse(formula)
                else:
                    formula = None

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            kegg = None if kegg.strip() == '' else kegg
            cas = None if cas.strip() == '' else cas

            yield CompoundEntry(id=compound_id, name=name,
                                formula=formula,
                                formula_neutral=formula_neutral,
                                charge=charge, kegg=kegg, cas=cas)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('Table 2')
        for i in range(1, sheet.nrows):
            (reaction_id, name, equation, _, genes, _, subsystem, ec,
                reversible) = sheet.row_values(i, end_colx=9)

            if reaction_id.strip() == '':
                continue

            if genes != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'[sb]\d+', genes))
            else:
                genes = None

            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = parse_sudensimple_reaction(equation)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            ec = None if ec.strip() == '' else ec

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem,
                                ec=ec)


class EColiTextbookImport(Importer):
    name = 'EColi_textbook'
    title = ('Escerichia coli Textbook (core) model (Excel format),'
             ' Orth et al., 2010')

    filename = 'ecoli_core_model.xls'

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('metabolites')
        for i in range(1, sheet.nrows):
            (compound_id, name, formula, charge, cas, formula_neutral,
                alt_names, kegg) = sheet.row_values(i, end_colx=8)

            if compound_id.strip() == '':
                continue

            # Skip compartmentalized compounds
            m = re.match(r'^(.*)\[.\]$', compound_id)
            if m:
                compound_id = m.group(1)

            name = None if name.strip() == '' else name
            formula = None if formula.strip() == '' else formula

            if formula_neutral.strip() != '':
                formula_neutral = Formula.parse(formula_neutral)
                formula = formula_neutral
            else:
                formula_neutral = None
                if formula.strip() != '':
                    formula = Formula.parse(formula)
                else:
                    formula = None

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            cas = None if cas.strip() == '' or cas == 'None' else cas
            kegg = None if kegg.strip() == '' else kegg

            yield CompoundEntry(id=compound_id, name=name,
                                formula=formula,
                                formula_neutral=formula_neutral,
                                charge=charge, kegg=kegg, cas=cas)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('reactions')
        for i in range(1, sheet.nrows):
            (reaction_id, name, equation, subsystem, ec, _, _, _, _, _,
                genes) = sheet.row_values(i, end_colx=11)

            if reaction_id.strip() == '':
                continue

            if genes.strip() != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'[sb]\d+', genes))
            else:
                genes = None

            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = parse_metnet_reaction(equation)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            ec = None if ec.strip() == '' else ec

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem, ec=ec)


class ImportSTMv1_0(Importer):
    name = 'STM_v1.0'
    title = ('Salmonella enterica STM_v1.0 (Excel format),'
             ' Thiele et al., 2011')

    filename = '1752-0509-5-8-s1.xlsx'

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())
        model.biomass_reaction = 'biomass_iRR1083_metals'

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('SI Tables - S2b - Metabolites')
        for i in range(2, sheet.nrows):
            _, compound_id, name, formula, charge, _, kegg, pubchem, chebi = (
                sheet.row_values(i, end_colx=9))

            if compound_id.strip() == '':
                continue

            name = None if name.strip() == '' else name.strip()

            if formula.strip() != '':
                formula = Formula.parse(formula)
            else:
                formula = None

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            kegg = None if kegg.strip() == '' else kegg

            yield CompoundEntry(id=compound_id, name=name,
                                formula=formula, charge=charge, kegg=kegg)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('SI Tables - S2a - Reactions')
        for i in range(4, sheet.nrows):
            reaction_id, name, equation, genes, _, subsystem = (
                sheet.row_values(i, end_colx=6))

            if reaction_id.strip() == '':
                continue

            name = None if name.strip() == '' else name.strip()

            if equation.strip() != '':
                equation = parse_sudensimple_reaction(equation,
                                                      arrow_irrev='-->')
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem

            if genes.strip() != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'STM\d+', genes))
            else:
                genes = None

            yield ReactionEntry(id=reaction_id, name=name,
                                genes=genes, equation=equation,
                                subsystem=subsystem)


class ImportiJN746(Importer):
    name = 'iJN746'
    title = ('Pseudomonas putida iJN746 (Excel format),'
             ' Nogales et al., 2011')

    filenames = ('1752-0509-2-79-s8.xls',
                 '1752-0509-2-79-s9.xls')

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:')
        for filename in self.filenames:
            print('- {}'.format(filename))

    def import_model(self, source):
        if os.path.isdir(source):
            compound_source = os.path.join(source, self.filenames[0])
            reaction_source = os.path.join(source, self.filenames[1])
        else:
            raise ParseError('Source must be a directory')

        self._compound_book = xlrd.open_workbook(compound_source)
        self._reaction_book = xlrd.open_workbook(reaction_source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._compound_book.sheet_by_name('Additional file 8')
        for i in range(1, sheet.nrows):
            (compound_id, name, formula, charge, cas, formula_neutral, _,
                kegg) = sheet.row_values(i, end_colx=8)

            if compound_id.strip() == '':
                continue

            name = None if name.strip() == '' else name

            if formula_neutral.strip() != '':
                formula_neutral = Formula.parse(formula_neutral)
                formula = formula_neutral
            else:
                formula_neutral = None
                if formula.strip() != '':
                    formula = Formula.parse(formula)
                else:
                    formula = None

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            kegg = None if kegg.strip() == '' else kegg

            if isinstance(cas, basestring) and cas.strip() in ('', 'None'):
                cas = None
            else:
                cas = str(cas)

            yield CompoundEntry(id=compound_id, name=name, formula=formula,
                                formula_neutral=formula_neutral,
                                charge=charge, kegg=kegg, cas=cas)

    def _read_reactions(self):
        sheet = self._reaction_book.sheet_by_name('Additional file 9')
        for i in range(1, sheet.nrows):
            reaction_id, name, equation, subsystem, ec, _, genes = (
                sheet.row_values(i, end_colx=7))

            if reaction_id.strip() == '':
                continue

            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = parse_metnet_reaction(equation)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            ec = None if ec.strip() == '' else ec

            if genes.strip() != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'PP_\d+', genes))
            else:
                genes = None

            yield ReactionEntry(id=reaction_id, name=name,
                                genes=genes, equation=equation,
                                subsystem=subsystem, ec=ec)


class ImportiJP815(Importer):
    name = 'iJP815'
    title = ('Pseudomonas putida iJP815 (Excel format),'
             ' Puchalka et al., 2008')

    filename = 'journal.pcbi.1000210.s011.XLS'

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('Metabolites')
        for i in range(1, sheet.nrows):
            compound_id, name = sheet.row_values(i, end_colx=2)

            if compound_id.strip() == '':
                continue

            # KEGG ID id encoded in the compound ID
            compartment, compound_id = (
                re.match(r'^(E|I)(C\d+)$', compound_id).groups())
            kegg = compound_id

            m = re.match('^(.*)\[.\]$', name)
            if m:
                name = m.group(1)

            yield CompoundEntry(id=compound_id, name=name, kegg=kegg)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('Reactions')
        for i in range(1, sheet.nrows):
            (reaction_id, name, equation, _, _, _, _, _, _, subsystem,
                genes) = sheet.row_values(i, end_colx=11)

            if reaction_id.strip() == '':
                continue

            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = parse_sudensimple_reaction(equation, '<==>', '-->')

                # Rebuild reaction with compartment information
                def translate(c, v):
                    compartment = 'e' if c.name[0] == 'E' else None
                    return Compound(c.name[1:], compartment=compartment), v

                left = (translate(c, v) for c, v in equation.left)
                right = (translate(c, v) for c, v in equation.right)
                equation = Reaction(equation.direction, left, right)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem

            if genes.strip() != '':
                genes = frozenset('PP_'+m.group(1)
                                  for m in re.finditer(r'PP(\d+)', genes))
            else:
                genes = None

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem)


class ImportiSyn731(Importer):
    name = 'iSyn731'
    title = ('Synechocystis sp. PCC 6803 iSyn731 (Excel format),'
             ' Saha et al., 2012')

    filename = 'journal.pone.0048285.s001.XLSX'

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())
        model.biomass_reaction = 'Biomass_Hetero'

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('Metabolites')
        for i in range(1, sheet.nrows):
            compound_id, name, formula, charge, kegg = (
                sheet.row_values(i))

            if compound_id.strip() == '':
                continue

            # Fixup model errors
            #m = re.match(r'^(.*)(C\d{5})$', formula_neutral)
            #if m:
            #    formula_neutral = m.group(1)
            #    kegg = m.group(2)

            name = None if name == '' else name
            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            #if formula_neutral.strip() != '':
            #    formula_neutral = Formula.parse(formula_neutral)
            #    formula = formula_neutral
            #else:
            #    formula_neutral = None
            if formula.strip() != '':
                formula = re.sub(r'-', 'noformula', formula)
                formula = re.sub(r'noformula', ' ', formula)
                m = re.match(r'^(.*)-\d$', formula)
                if m:
                    formula = m.group(1)
                formula = Formula.parse(formula)
            else:
                formula = None

            if kegg != 0 and kegg.strip() != '':
                # Discard secondary KEGG IDs
                kegg = kegg.split('|')[0]
            else:
                kegg = None

            yield CompoundEntry(id=compound_id, name=name, formula=formula,
                                charge=charge, kegg=kegg)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('Model')
        for i in range(2, sheet.nrows):
            reaction_id, name, ec, genes, _, equation, subsystem = (
                sheet.row_values(i, end_colx=7))
            if reaction_id.strip() == 'EX_Arsenic acid':
                reaction_id = 'EX_Arsenic_acid'
            if reaction_id.strip() == '':
                continue

            name = None if name == '' else name

            if equation.strip() != '':
                # check that this works correctly. should substitute the => for a
                # space then =>. the double spaces should be ignored though.
                equation = re.sub(r'=>', ' =>', equation)
                equation = re.sub(r'< =>', '<=>', equation)
                equation = re.sub(r'\s+', ' ', equation)
                equation = re.sub(r'\(1\)\|', '(1) |', equation)
                equation = re.sub(r'\+\(1\)', '+ (1)', equation)
                equation = re.sub(r'\|\|', '|', equation)
                equation = modelseed.parse_reaction(equation)
            else:
                equation = None

            if genes != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'\w+\d+', genes))
            else:
                genes = None

            ec = ec if ec.strip() != '' and ec != 'Undetermined' else None

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem, ec=ec)


class ImportiCce806(Importer):
    name = 'iCce806'
    title = ('Cyanothece sp. ATCC 51142 iCce806 (Excel format),'
             ' Vu et al., 2012')

    filenames = ('journal.pcbi.1002460.s005.XLSX',
                 'journal.pcbi.1002460.s006.XLSX')

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:')
        for filename in self.filenames:
            print('- {}'.format(filename))

    def import_model(self, source):
        if os.path.isdir(source):
            reaction_source = os.path.join(source, self.filenames[0])
            compound_source = os.path.join(source, self.filenames[1])
        else:
            raise ParseError('Source must be a directory')

        self._compound_book = xlrd.open_workbook(compound_source)
        self._reaction_book = xlrd.open_workbook(reaction_source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())
        model.biomass_reaction = 'CyanoBM (average)'

        return model

    def _read_compounds(self):
        sheet = self._compound_book.sheet_by_name('Table S2')
        for i in range(1, sheet.nrows):
            compound_id, name, formula, charge, cas, formula_neutral, _, kegg = (
                sheet.row_values(i))

            if compound_id.strip() == '':
                continue

            # Fixup model errors
            m = re.match(r'^(.*)(C\d{5})$', formula_neutral)
            if m:
                formula_neutral = m.group(1)
                kegg = m.group(2)

            name = None if name == '' else name

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            #if formula_neutral.strip() != '':
            #    formula_neutral = formula.parse(formula_neutral)
            #    #formula = formula_neutral
            #else:
            #    formula_neutral = None
           # print(formula_neutral)
            if formula.strip() != '':
                m = re.match(r'^(.*)-\d$', formula)
                if m:
                    formula = m.group(1)
                formula = Formula.parse(formula)
            else:
                formula = None

            kegg = None if kegg == '' else kegg

            if cas.strip() != '' and cas.strip() != 'None':
                cas = cas.strip()
            else:
                cas = None

            yield CompoundEntry(id=compound_id, name=name, formula=formula,
                                formula_neutral=formula, charge=charge,
                                kegg=kegg, cas=cas)

    def _read_reactions(self):
        sheet = self._reaction_book.sheet_by_name('S1 - Reactions')
        for i in range(1, sheet.nrows):
            reaction_id, name, equation, _, genes, subsystem, ec = (
                sheet.row_values(i, end_colx=7))

            if reaction_id.strip() == '':
                continue
            if reaction_id.strip() == 'Notes:':
                continue
            if reaction_id.strip() == 'Abbreviation':
                continue
            if reaction_id.strip() == 'AL':
                continue
            if reaction_id.strip() == 'LL':
                continue
            if reaction_id.strip() == 'Column headings':
                continue
            if reaction_id.strip() == 'Column H through K':
                continue
            if reaction_id.strip() == 'Column H':
                continue
            if reaction_id.strip() == 'Column I':
                continue
            if reaction_id.strip() == 'Column J':
                continue
            if reaction_id.strip() == 'Column K':
                continue
            # Fixup model errors

            name = None if name == '' else name

            if equation.strip() != '':
                equation = re.sub(r'\s+', ' ', equation)
                equation = parse_metnet_reaction(equation)
            else:
                equation = None
            if genes != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'cce\_[0-9]+', genes))
            else:
                genes = None

            m = re.match(r'EC-(.*)', ec)
            if m:
                ec = m.group(1)
                if ec == 'Undetermined':
                    ec = None
            else:
                ec = None

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem, ec=ec)


class ImportGSMN_TB(Importer):
    name = 'GSMN-TB'
    title = ('Mycobacterium tuberculosis GSMN-TB (Excel format),'
             ' Beste et al., 2007')

    filenames = ('gb-2007-8-5-r89-s4.xls',
                 'gb-2007-8-5-r89-s6.xls')

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:')
        for filename in self.filenames:
            print('- {}'.format(filename))

    def import_model(self, source):
        if os.path.isdir(source):
            reaction_source = os.path.join(source, self.filenames[0])
            compound_source = os.path.join(source, self.filenames[1])
        else:
            raise ParseError('Source must be a directory')

        self._compound_book = xlrd.open_workbook(compound_source)
        self._reaction_book = xlrd.open_workbook(reaction_source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._compound_book.sheet_by_name('File 6')
        for i in range(2, sheet.nrows):
            compound_id, name = sheet.row_values(i, end_colx=2)

            if compound_id.strip() == '':
                continue

            # Skip compartmentalized compounds
            m = re.match(r'^(.*)\[.\]$', compound_id)
            if m:
                compound_id = m.group(1)

            name = None if name.strip() == '' else name
            #formula = None if formula.strip() == '' else formula

            #if formula_neutral.strip() != '':
            #    formula_neutral = Formula.parse(formula_neutral)
           #     formula = formula_neutral
           # else:
            #    formula_neutral = None
            #    if formula.strip() != '':
            #        formula = Formula.parse(formula)
            #    else:
            #        formula = None

            #try:
            #    charge = None if charge == '' else int(charge)
            #except ValueError:
            #    charge = None

           # cas = None if cas.strip() == '' or cas == 'None' else cas
           # kegg = None if kegg.strip() == '' else kegg

            yield CompoundEntry(id=compound_id, name=name)

        def create_missing(compound_id, name=None):
            if name is None:
                name = compound_id
            return CompoundEntry(id=compound_id, name=name)

        # Generate missing compounds
        yield create_missing('MBT-HOLO', 'Mycobactin-Holo')
        yield create_missing('TAGbio', 'Triacylglycerol-bio')
        yield create_missing('TREHALOSEMONOMYCOLATE(CY)',
                             'Trehalosemonomycolate(cy)')
        yield create_missing('MYCOTHIOL-S-CONJUGATE', 'Mycothiol-S-conjugate')
        yield create_missing('N-ACETYL-S-CONJUGATE', 'N-actyl-s-conjugate')
        yield create_missing('(n-1)POLYP', '(n-1)Polyphosphate')
        yield create_missing('(n)POLYP', '(n)Polyphosphate')
        yield create_missing('HYDROXYGLU', 'Hydroxy-glutamate')
        yield create_missing('OCTANOYL-ACP', 'OCTANOYL-Acyl-Carrier-Protein')
        yield create_missing('UDP[NAM:NGM]ALA')
        yield create_missing('UDP[NAM:NGM]ALAGLU')
        yield create_missing('UDP[NAM:NGM]AGMDAPIMAA')
        yield create_missing('UDP[NAM:NGM]AGMDAPIM')
        yield create_missing('MOLYBDENUM-CO', 'MOLYBDENUM-COfactor')
        yield create_missing('CISACONITATE')
        yield create_missing('ELECTROPHILE-X')
        yield create_missing('H2X')
        yield create_missing('TAGcat', 'TracylGlycerol Cat')
        yield create_missing('METHYLISOCITRATE', 'Methyl-Isocitrate')
        yield create_missing('ARA[1]GALACTANDPP')
        yield create_missing('MPM', 'Methyl Pyrimidine')
        yield create_missing('APO-LIPO', 'apo-lipo-amide')
        yield create_missing('BIOMASSxt', 'Biomass extracellular')
        yield create_missing('GLUCAN', 'Glucanate')

    def _read_reactions(self):
        sheet = self._reaction_book.sheet_by_name('File 4')
        for i in range(5, sheet.nrows):
            (reaction_id, equation, fluxbound, _, ec, genes, name,
                subsystem) = sheet.row_values(i, end_colx=8)


            if reaction_id.startswith('%') or reaction_id.strip() == '':
                continue
            #if not reaction_id.stripstartswith('%'):
            #    continue
            #print(equation)
            if genes.strip() != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'Rv\w+', genes))
            else:
                genes = None

            name = None if name.strip() == '' else name
            if equation.strip() != '':
                equation = re.sub(r'\s+', ' ', equation)
                equation = parse_sudensimple_reaction(equation, '=')
                rdir = Reaction.Bidir if fluxbound != 0 else Reaction.Right
                equation = Reaction(rdir, equation.left, equation.right)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            ec = None if ec.strip() == '' else ec

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem, ec=ec)


class ImportiNJ661(Importer):
    name = 'iNJ661'
    title = ('Mycobacterium tuberculosis iNJ661 (Excel format),'
             ' Jamshidi et al., 2007')

    filename = '1752-0509-1-26-s5.xls'

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            self.title, self._read_compounds(), self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('metabolites')
        for i in range(1, sheet.nrows):
            compound_id, name, formula, charge = sheet.row_values(
                i, end_colx=4)

            if compound_id.strip() == '':
                continue

            name = name if name.strip() != '' else None

            if formula.strip() != '':
                formula = Formula.parse(formula)
            else:
                formula = None

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            yield CompoundEntry(id=compound_id, name=name, formula=formula,
                                charge=charge)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('iNJ661')
        for i in range(5, sheet.nrows):
            reaction_id, name, equation, _, subsystem, _, genes = (
                sheet.row_values(i, end_colx=7))

            if reaction_id.strip() == '':
                continue

            if genes.strip() != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'Rv\w+', genes))
            else:
                genes = None

            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = parse_metnet_reaction(equation)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem)


class ImportGenericiNJ661mv(Importer):
    """Importer for the model iNJ661m and iNJ661v

    For models of Mycobacterium tuberculosis iNJ661m/v (Excel format),
    Fang et al., 2010.
    """

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model_named(self, name, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)

        model = MetabolicModel(
            name, self._read_compounds(), self._read_reactions())
        model.biomass_reaction = 'biomass_Mtb_9_60atp_test_NOF'

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('metabolites')
        for i in range(1, sheet.nrows):
            compound_id, name, formula = sheet.row_values(i, end_colx=3)

            # Skip compartmentalized compounds
            compound_id = re.match(r'^(.*)\[.\]$', compound_id).group(1)
            name = name if name.strip() != '' else None

            if formula.strip() != '':
                formula = Formula.parse(formula)
            else:
                formula = None

            yield CompoundEntry(id=compound_id, name=name, formula=formula)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('reactions')
        for i in range(1, sheet.nrows):
            reaction_id, name, equation, genes, _, subsystem = (
                sheet.row_values(i, end_colx=6))

            if reaction_id.strip() == '':
                continue

            if genes.strip() != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'Rv\w+', genes))
            else:
                genes = None

            name = None if name.strip() == '' else name.strip()

            # Biomass reaction is not specified in this table
            if (equation.strip() != '' and
                    reaction_id != 'biomass_Mtb_9_60atp_test_NOF'):
                # Remove multiple adjacent whitespace characters
                equation = re.sub(r'\s+', ' ', equation)
                equation = parse_sudensimple_reaction(
                    equation, arrow_rev='<=>', arrow_irrev='->')
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem)


class ImportiNJ661m(ImportGenericiNJ661mv):
    name = 'inj661m'
    title = ('Mycobacterium tuberculosis iNJ661m (Excel format),'
        ' Fang et al., 2010')
    filename = '1752-0509-4-160-s3.xls'

    def import_model(self, source):
        return self.import_model_named(self.title, source)


class ImportiNJ661v(ImportGenericiNJ661mv):
    name = 'inj661v'
    title = ('Mycobacterium tuberculosis iNJ661v (Excel format),'
        ' Fang et al., 2010')
    filename = '1752-0509-4-160-s5.xls'

    def import_model(self, source):
        return self.import_model_named(self.title, source)


class ImportShewanellaOng(Importer):
    """Generic importer for four models published in Ong et al., 2014

    Generic importer for the models iMR1_799, iMR4_812, iW3181_789 and
    iOS217_672 from Ong et al. 2014. "Comparisons of Shewanella Strains Based
    on Genome Annotations, Modeling, and Experiments." BMC Systems Biology 8
    (1). BMC Systems Biology: 1-11. doi:10.1186/1752-0509-8-31.
    """

    filename = '1752-0509-8-31-s2.xlsx'
    biomass_names = (
        'SO_BIOMASSMACRO_DM_NOATP2',
        'MR4_BIOMASSMACRO_DM_NOATP2',
        'W3181_BIOMASSMACRO_DM_NOATP2',
        'Sden_BIOMASSMACRO_DM_NOATP2',
        'Core_BIOMASSMACRO_DM_NOATP2'
    )

    def help(self):
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model_named(self, name, col_index, source):
        if os.path.isdir(source):
            source = os.path.join(source, self.filename)

        self._book = xlrd.open_workbook(source)
        self._col_index = col_index

        model = MetabolicModel(
            name, self._read_compounds(), self._read_reactions())
        model.biomass_reaction = self.biomass_names[col_index]

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('S3-Metabolites')
        for i in range(1, sheet.nrows):
            (compound_id, _, _, _, _, _, name, formula_neutral, formula,
                charge, _, kegg, cas) = sheet.row_values(i, end_colx=13)

            # Remove compartmentalization of compounds
            compound_id = re.match(r'^(.*)\[.\]$', compound_id).group(1)
            name = name if name.strip() != '' else None

            if formula_neutral.strip() != '':
                formula_neutral = Formula.parse(formula_neutral)
                formula = formula_neutral
            else:
                formula_neutral = None
                if formula.strip() != '':
                    formula = Formula.parse(formula)
                else:
                    formula = None

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            kegg = None if kegg.strip() == '' else kegg

            if isinstance(cas, basestring) and cas.strip() in ('', 'None'):
                cas = None
            else:
                cas = str(cas)

            yield CompoundEntry(id=compound_id, name=name, formula=formula,
                                formula_neutral=formula_neutral, charge=charge,
                                kegg=kegg, cas=cas)

    def _read_reactions(self):
        sheet = self._book.sheet_by_name('S2-Reactions')
        for i in range(2, sheet.nrows):
            reaction_id, _, _, _, _, _, name, equation = sheet.row_values(
                i, end_colx=8)

            if reaction_id.strip() == '':
                continue

            # Whether the reaction is present in this model
            model_presence = bool(sheet.row_values(
                i, start_colx=1, end_colx=6)[self._col_index])

            if not model_presence:
                # TODO load the complete reaction list and use the model subset
                # feature of the native format to restrict the model.
                continue

            # Genes
            model_genes = sheet.row_values(
                i, start_colx=13, end_colx=18)[self._col_index].strip()

            # Fixup compound names in reactions
            def translate(s):
                s = s.lower()
                m = re.match(r'^(.*)_e$', s)
                if m:
                    s = m.group(1)

                if s == 'aaacoa':
                    s = 'aacoa'

                if s in ('fdxr-4:2', 'fdxo-4:2'):
                    s = s.replace(':', '_')

                if s in ('q8', 'q8h2'):
                    s = 'ub' + s
                return s

            # Reaction equation
            if equation.strip() != '':
                equation = parse_metnet_reaction(equation)
                equation = equation.translated_compounds(translate)
            else:
                equation = None

            if model_genes != '':
                genes = frozenset(model_genes.split())
            else:
                genes = None

            subsystem = sheet.cell_value(i, 18)

            name = None if name.strip() == '' else name.strip()
            subsystem = None if subsystem.strip() == '' else subsystem

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, subsystem=subsystem)


class ImportiMR1_799(ImportShewanellaOng):
    name = 'imr1_799'
    title = ('Shewanella oneidensis MR-1 iMR1_799 (Excel format),'
        ' Ong et al., 2014')

    def import_model(self, source):
        return self.import_model_named(self.title, 0, source)


class ImportiMR4_812(ImportShewanellaOng):
    name = 'imr4_812'
    title = ('Shewanella sp. MR-4 iMR4_812 (Excel format),'
        ' Ong et al., 2014')

    def import_model(self, source):
        return self.import_model_named(self.title, 1, source)


class ImportiW3181_789(ImportShewanellaOng):
    name = 'iw3181_789'
    title = ('Shewanella sp. W3-18-1 iW3181_789 (Excel format),'
        ' Ong et al., 2014')

    def import_model(self, source):
        return self.import_model_named(self.title, 2, source)


class ImportiOS217_672(ImportShewanellaOng):
    name = 'ios217_672'
    title = ('Shewanella denitrificans OS217 iOS217_672 (Excel format),'
        ' Ong et al., 2014')

    def import_model(self, source):
        return self.import_model_named(self.title, 3, source)


class ImportModelSEED(Importer):
    """Read metabolic model for a ModelSEED model"""

    name = 'ModelSEED'
    title = 'ModelSEED model (Excel format)'

    def help(self):
        print('Source must contain the model definition in Excel format\n'
              ' and a PTT file for mapping PEG gene names.'
              'Expected files in source directory:\n'
              '- Seed*.xls\n'
              '- NC_*.ptt')

    def import_model(self, source):
        if not os.path.isdir(source):
            raise ParseError('Source must be a directory')

        excel_sources = glob.glob(os.path.join(source, 'Seed*.xls'))
        if len(excel_sources) == 0:
            raise ParseError('No .xls file found in source directory')
        elif len(excel_sources) > 1:
            raise ParseError(
                'More than one .xls file found in source directory')

        ptt_sources = glob.glob(os.path.join(source, '*.ptt'))
        if len(ptt_sources) == 0:
            raise ParseError('No .ptt file found in source directory')
        elif len(ptt_sources) > 1:
            raise ParseError(
                'More than one .ptt file found in source directory')

        self._book = xlrd.open_workbook(excel_sources[0])

        with open(ptt_sources[0], 'r') as ptt_file:
            # Read mapping from location to gene ID from PTT file
            location_mapping = {}
            for i in range(3):
                ptt_file.readline() # Skip headers
            for row in csv.reader(ptt_file, delimiter='\t'):
                location, direction, _, _, _, gene_id = row[:6]

                start, stop = re.match(r'^(\d+)\.\.(\d+)$', location).groups()
                start, stop = int(start), int(stop)

                location_mapping[start, stop, direction] = gene_id

        # Read mapping from PEG to gene ID
        peg_mapping = {}
        sheet = self._book.sheet_by_name('Genes')
        for i in range(1, sheet.nrows):
            gene_id, gene_type, _, start, stop, direction = sheet.row_values(
                i, end_colx=6)

            m = re.match(r'^.*(peg.\d+)$', gene_id)
            if not m:
                continue

            peg_id = m.group(1)
            direction = '+' if direction == 'for' else '-'
            start, stop = int(start), int(stop)

            peg_mapping[peg_id] = location_mapping[start, stop, direction]

        model = MetabolicModel(
            'ModelSEED model', self._read_compounds(),
            self._read_reactions(peg_mapping))

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('Compounds')
        for i in range(1, sheet.nrows):
            compound_id, name, alt_name, formula, charge, _ = (
                sheet.row_values(i, end_colx=6))

            if compound_id.strip() == '':
                continue

            charge = int(charge) if charge != '' else None

            if formula.strip() != '':
                formula = formula.strip()
            else:
                formula = None

            yield CompoundEntry(id=compound_id, name=name, formula=formula,
                                charge=charge)

    def _read_reactions(self, peg_mapping):
        sheet = self._book.sheet_by_name('Reactions')
        for i in range(1, sheet.nrows):
            reaction_id, name, equation, _, ec_list, _, _, pegs = (
                sheet.row_values(i, end_colx=8))

            if reaction_id.strip() == '':
                continue

            name = name if name.strip() != '' else None

            if equation != '' and 'NONE' not in equation:
                equation = modelseed.parse_reaction(equation)
            else:
                continue

            if pegs != '':
                pegs = frozenset(m.group(0) for m in
                                  re.finditer(r'peg\.(\d+)', pegs))
                genes = frozenset(peg_mapping[p] for p in pegs)
            else:
                pegs = None
                genes = None

            if ec_list == '':
                ec_list = None
                ec = None
            elif '|' in ec_list:
                ec_list = frozenset(ec.strip() for ec in ec_list.split('|')
                                    if ec.strip() != '')
                ec = next(iter(ec_list))
            else:
                ec_list = frozenset([ec_list])
                ec = next(iter(ec_list))

            yield ReactionEntry(id=reaction_id, name=name, genes=genes,
                                equation=equation, ec=ec)


class SBMLImporter(Importer):
    """Base importer for reading metabolic model from an SBML file"""

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


class SBMLStrictImporter(SBMLImporter):
    """Read metabolic model from an SBML file using strict parser"""

    name = 'SBML-strict'
    title = 'SBML model (strict)'

    def _open_reader(self, f):
        return sbml.SBMLReader(f, strict=True, ignore_boundary=True)


class SBMLNonstrictImporter(SBMLImporter):
    """Read metabolic model from an SBML file using non-strict parser"""

    name = 'SBML'
    title = 'SBML model (non-strict)'

    def _open_reader(self, f):
        return sbml.SBMLReader(f, strict=False, ignore_boundary=True)

    def import_model(self, source):
        model = super(SBMLNonstrictImporter, self).import_model(source)

        biomass_reaction = None
        objective_reactions = set()
        for reaction in self._reader.reactions:
            # Check whether species multiple times
            compounds = set()
            for c, v in reaction.equation.left:
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

        if len(objective_reactions) == 1:
            biomass_reaction = next(iter(objective_reactions))
            logger.info('Detected biomass reaction: {}'.format(
                biomass_reaction))

        model = MetabolicModel(
            model.name,
            self._convert_compounds(model.compounds.itervalues()),
            self._convert_reactions(model.reactions.itervalues()))
        model.biomass_reaction = biomass_reaction

        return model

    def _convert_compounds(self, compounds):
        """Convert SBML species entries to compounds"""
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

    def _convert_reactions(self, reactions):
        """Convert SBML reaction entries to reactions"""
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

            yield ReactionEntry(**properties)
