
"""Functions related to loading models"""

from __future__ import print_function

import os
import re
import csv
from collections import namedtuple

import xlrd
from metnet.datasource.misc import (parse_metnet_reaction,
                                    parse_sudensimple_reaction)
from metnet.datasource import modelseed
from metnet.formula import Formula
from metnet.reaction import Reaction, Compound


CompoundEntry = namedtuple('CompoundEntry',
                           ['id', 'name', 'formula', 'formula_neutral',
                            'charge', 'kegg', 'cas'])
ReactionEntry = namedtuple('ReactionEntry',
                           ['id', 'name', 'genes', 'equation',
                            'subsystem', 'ec'])


class Importer(object):
    """Base importer class"""


class MetabolicModel(object):
    def __init__(self, name, compounds, reactions):
        self._name = name
        self._compounds = dict(compounds)
        self._reactions = dict(reactions)

        self._genes = set()
        for r in self._reactions.itervalues():
            if r.genes is not None:
                self._genes.update(r.genes)

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

    def print_summary(self):
        """Print model summary"""
        print('Model: {}'.format(self.name))
        print('- Compounds: {}'.format(len(self.compounds)))
        print('- Reactions: {}'.format(len(self.reactions)))
        print('- Genes: {}'.format(len(self.genes)))

    def check_reaction_compounds(self):
        """Check that reaction compounds are defined in the model"""
        for reaction in self.reactions.itervalues():
            if reaction.equation is not None:
                for compound, value in reaction.equation.compounds:
                    if compound.name not in self.compounds:
                        return reaction.id, compound.name
        return None, None


class ParseError(Exception):
    """Exception used to signal a parsing error"""


class ImportiMA945(Importer):
    name = 'ima945'
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
            'iMA945', self._read_compounds(), self._read_reactions())

        reaction_id, compound_name = model.check_reaction_compounds()
        if compound_name is not None:
            raise ParseError(
                'Compound {}, {} not defined in compound table'.format(
                    reaction_id, compound_name))
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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=formula_neutral,
                                  charge=charge, kegg=kegg, cas=None)
            yield compound_id, entry

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

            yield reaction_id, ReactionEntry(id=reaction_id, name=name,
                                             genes=genes, equation=equation,
                                             subsystem=None, ec=None)


class ImportiRR1083(Importer):
    name = 'irr1083'
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
            'iRR1083', self._read_compounds(), self._read_reactions())

        reaction_id, compound_name = model.check_reaction_compounds()
        if compound_name is not None:
            raise ParseError(
                'Compound {}, {} not defined in compound table'.format(
                    reaction_id, compound_name))

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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula_neutral,
                                  formula_neutral=formula_neutral,
                                  charge=charge, kegg=kegg, cas=None)
            yield compound_id, entry

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

            entry = ReactionEntry(id=reaction_id, name=name, genes=genes,
                                  equation=equation, subsystem=None, ec=None)
            yield reaction_id, entry


class ImportiJO1366(Importer):
    name = 'ijo1366'
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
            'iJO1366', self._read_compounds(), self._read_reactions())

        reaction_id, compound_name = model.check_reaction_compounds()
        if compound_name is not None:
            raise ParseError(
                'Compound {}, {} not defined in compound table'.format(
                    reaction_id, compound_name))
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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=formula_neutral,
                                  charge=charge, kegg=kegg, cas=cas)
            yield compound_id, entry

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

            entry = ReactionEntry(id=reaction_id, name=name, genes=genes,
                                  equation=equation, subsystem=subsystem,
                                  ec=ec)
            yield reaction_id, entry


class EColiTextbookImport(Importer):
    name = 'ecoli_textbook'
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
            'EColi_textbook', self._read_compounds(), self._read_reactions())

        reaction_id, compound_name = model.check_reaction_compounds()
        if compound_name is not None:
            raise ParseError('Compound {}, {} not defined in compound table'.format(
                reaction_id, compound_name))

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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=formula_neutral,
                                  charge=charge, kegg=kegg, cas=cas)
            yield compound_id, entry

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

            entry = ReactionEntry(id=reaction_id, name=name, genes=genes,
                                  equation=equation, subsystem=subsystem,
                                  ec=ec)
            yield reaction_id, entry


class ImportSTMv1_0(Importer):
    name = 'stm_v1.0'
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
            'STM_v1.0', self._read_compounds(), self._read_reactions())

        reaction_id, compound_name = model.check_reaction_compounds()
        if compound_name is not None:
            raise ParseError(
                'Compound {}, {} not defined in compound table'.format(
                    reaction_id, compound_name))
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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=None,
                                  charge=charge, kegg=kegg, cas=None)
            yield compound_id, entry

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

            entry = ReactionEntry(id=reaction_id, name=name,
                                  genes=genes, equation=equation,
                                  subsystem=subsystem, ec=None)
            yield reaction_id, entry


class ImportiJN746(Importer):
    name = 'ijn746'
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
            'iJN746', self._read_compounds(), self._read_reactions())

        reaction_id, compound_name = model.check_reaction_compounds()
        if compound_name is not None:
            raise ParseError(
                'Compound {}, {} not defined in compound table'.format(
                    reaction_id, compound_name))
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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=formula_neutral,
                                  charge=charge, kegg=kegg, cas=cas)
            yield compound_id, entry

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

            entry = ReactionEntry(id=reaction_id, name=name,
                                  genes=genes, equation=equation,
                                  subsystem=subsystem, ec=ec)
            yield reaction_id, entry



class ImportiJP815(Importer):
    name = 'ijp815'
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
            'iJP815', self._read_compounds(), self._read_reactions())

        reaction_id, compound_name = model.check_reaction_compounds()
        if compound_name is not None:
            raise ParseError('Compound {}, {} not defined in compound table'.format(
                reaction_id, compound_name))

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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=None,
                                  formula_neutral=None,
                                  charge=None, kegg=kegg, cas=None)
            yield compound_id, entry

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

            entry = ReactionEntry(id=reaction_id, name=name,
                                  genes=genes, equation=equation,
                                  subsystem=subsystem, ec=None)
            yield reaction_id, entry


class ImportiSyn731(Importer):
    name = 'isyn731'
    title = ('Synechocystis sp. PCC 6803 iSyn731 (Excel format),'
             ' Saha et al., 2008')

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
            'iSyn731', self._read_compounds(), self._read_reactions())

        reaction_id, compound_name = model.check_reaction_compounds()
        if compound_name is not None:
            raise ParseError(
                'Compound {}, {} not defined in compound table'.format(
                    reaction_id, compound_name))
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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=None,
                                  charge=charge, kegg=kegg, cas=None)
            yield compound_id, entry

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

            entry = ReactionEntry(id=reaction_id, name=name, genes=genes,
                                  equation=equation, subsystem=subsystem,
                                  ec=ec)
            yield reaction_id, entry


def read_cyanotheca_vu(reaction_excel_path, compound_excel_path):
    """Read compounds and reactions for cyanotheca Vu model"""

    reaction_book = xlrd.open_workbook(reaction_excel_path)
    compound_book = xlrd.open_workbook(compound_excel_path)

    def read_compounds(sheet):
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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=formula,
                                  charge=charge, kegg=kegg, cas=cas)
            yield compound_id, entry

    def read_reactions(sheet):
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

            yield reaction_id, ReactionEntry(id=reaction_id, name=name,
                                             genes=genes, equation=equation,
                                             subsystem=subsystem, ec=ec)
    model = MetabolicModel(
        'cyanotheca_vu',
        read_compounds(compound_book.sheet_by_name('Table S2')),
        read_reactions(reaction_book.sheet_by_name('S1 - Reactions')))

    reaction_id, compound_name = model.check_reaction_compounds()
    if compound_name is not None:
        raise ParseError('Compound {}, {} not defined in compound table'.format(
            reaction_id, compound_name))

    return model


def read_mtuberc_beste(reaction_excel_path, compound_excel_path):
    """Read compounds and reactions for M. tuberculosis Beste model"""

    reaction_book = xlrd.open_workbook(reaction_excel_path)
    compound_book = xlrd.open_workbook(compound_excel_path)

    def read_compounds(sheet):
        for i in range(2, sheet.nrows):
            (compound_id, name) = sheet.row_values(i, end_colx=2)

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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=None,
                                  formula_neutral=None,
                                  charge=None, kegg=None, cas=None)
            yield compound_id, entry

        # Generate missing compounds
        yield 'MBT-HOLO', CompoundEntry(id='MBT-HOLO', name='Mycobactin-Holo',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'TAGbio', CompoundEntry(id='TAGbio', name='Triacylglycerol-bio',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'TREHALOSEMONOMYCOLATE(CY)', CompoundEntry(id='TREHALOSEMONOMYCOLATE(CY)', name='Trehalosemonomycolate(cy)',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'MYCOTHIOL-S-CONJUGATE', CompoundEntry(id='MYCOTHIOL-S-CONJUGATE', name='Mycothiol-S-conjugate',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'N-ACETYL-S-CONJUGATE', CompoundEntry(id='N-ACETYL-S-CONJUGATE', name='N-actyl-s-conjugate',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield '(n-1)POLYP', CompoundEntry(id='(n-1)POLYP', name='(n-1)Polyphosphate',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield '(n)POLYP', CompoundEntry(id='(n)POLYP', name='(n)Polyphosphate',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'HYDROXYGLU', CompoundEntry(id='HYDROXYGLU', name='Hydroxy-glutamate',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'OCTANOYL-ACP', CompoundEntry(id='OCTANOYL-ACP', name='OCTANOYL-Acyl-Carrier-Protein',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'UDP[NAM:NGM]ALA', CompoundEntry(id='UDP[NAM:NGM]ALA', name='UDP[NAM:NGM]ALA',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'UDP[NAM:NGM]ALAGLU', CompoundEntry(id='UDP[NAM:NGM]ALAGLU', name='UDP[NAM:NGM]ALAGLU',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'UDP[NAM:NGM]AGMDAPIMAA', CompoundEntry(id='UDP[NAM:NGM]AGMDAPIMAA', name='UDP[NAM:NGM]AGMDAPIMAA',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'UDP[NAM:NGM]AGMDAPIM', CompoundEntry(id='UDP[NAM:NGM]AGMDAPIM', name='UDP[NAM:NGM]AGMDAPIM',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'MOLYBDENUM-CO', CompoundEntry(id='MOLYBDENUM-CO', name='MOLYBDENUM-COfactor',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'CISACONITATE', CompoundEntry(id='CISACONITATE', name='CISACONITATE',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'ELECTROPHILE-X', CompoundEntry(id='ELECTROPHILE-X', name='ELECTROPHILE-X',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'H2X', CompoundEntry(id='H2X', name='H2X',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'TAGcat', CompoundEntry(id='TAGcat', name='TracylGlycerol Cat',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'METHYLISOCITRATE', CompoundEntry(id='METHYLISOCITRATE', name='Methyl-Isocitrate',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'ARA[1]GALACTANDPP', CompoundEntry(id='ARA[1]GALACTANDPP', name='ARA[1]GALACTANDPP',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'MPM', CompoundEntry(id='MPM', name='Methyl Pyrimidine',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'APO-LIPO', CompoundEntry(id='APO-LIPO', name='apo-lipo-amide',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'BIOMASSxt', CompoundEntry(id='BIOMASSxt', name='Biomass extracellular',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)
        yield 'GLUCAN', CompoundEntry(id='GLUCAN', name='Glucanate',formula=None, formula_neutral=None, charge=None, kegg=None, cas=None)


    def read_reactions(sheet):
        for i in range(5, sheet.nrows):
            (reaction_id, equation, fluxbound, _, ec, genes, name, subsystem) = sheet.row_values(i, end_colx=8)


            if reaction_id.startswith('%') or reaction_id.strip() == '':
                continue
            #if not reaction_id.stripstartswith('%'):
            #    continue
            #print(equation)
            genes = None if genes.strip() == '' else frozenset(m.group(0) for m in re.finditer(r'[Rv]\d+[c]', genes))
            name = None if name.strip() == '' else name
            #equation = None if equation.strip() == '' else parse_sudensimple_reaction_fluxdef(equation, fluxbound)
            if equation.strip() != '':
                equation = re.sub(r'\s+', ' ', equation)
                equation = parse_sudensimple_reaction(equation, '=')
                direction = Reaction.Bidir if fluxbound != 0 else Reaction.Right
                equation = Reaction(direction, equation.left, equation.right)
            else:
                equation = None
            subsystem = None if subsystem.strip() == '' else subsystem
            yield reaction_id, ReactionEntry(id=reaction_id, name=name, genes=genes, equation=equation,
                                             subsystem=subsystem, ec=ec)

    model = MetabolicModel(
        'M_tuberculosis_beste',
        read_compounds(compound_book.sheet_by_name('File 6')),
        read_reactions(reaction_book.sheet_by_name('File 4')))

    reaction_id, compound_name = model.check_reaction_compounds()
    if compound_name is not None:
        raise ParseError('Compound {}, {} not defined in compound table'.format(
            reaction_id, compound_name))

    return model


def read_mtuberc_jamshidi(excel_path):
    """Read compounds and reactions for M. tuberculosis Jamshidi (iNJ661)"""

    book = xlrd.open_workbook(excel_path)

    def read_compounds(sheet):
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

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=None,
                                  charge=charge, kegg=None, cas=None)
            yield compound_id, entry

    def read_reactions(sheet):
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
            entry = ReactionEntry(id=reaction_id, name=name,
                                  genes=genes, equation=equation,
                                  subsystem=subsystem, ec=None)
            yield reaction_id, entry

    model = MetabolicModel(
        'M_tuberculosis_jamshidi',
        read_compounds(book.sheet_by_name('metabolites')),
        read_reactions(book.sheet_by_name('iNJ661')))

    reaction_id, compound_name = model.check_reaction_compounds()
    if compound_name is not None:
        raise ParseError('Compound {}, {} not defined in compound table'.format(
            reaction_id, compound_name))

    return model


def read_mtuberc_fang(name, excel_path):
    """Read compounds and reactions for M. tuberculosis FANG iNJ661m/v model"""

    book = xlrd.open_workbook(excel_path)

    def read_compounds(sheet):
        for i in range(1, sheet.nrows):
            compound_id, name, formula = sheet.row_values(i, end_colx=3)

            # Skip compartmentalized compounds
            compound_id = re.match(r'^(.*)\[.\]$', compound_id).group(1)
            name = name if name.strip() != '' else None

            if formula.strip() != '':
                formula = Formula.parse(formula)
            else:
                formula = None

            entry = CompoundEntry(id=compound_id, name=name,
                                  formula=formula,
                                  formula_neutral=None,
                                  charge=None, kegg=None, cas=None)
            yield compound_id, entry

    def read_reactions(sheet):
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

            name = None if name.strip() == '' else name

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

            entry = ReactionEntry(id=reaction_id, name=name,
                                  genes=genes, equation=equation,
                                  subsystem=subsystem, ec=None)
            yield reaction_id, entry

    model = MetabolicModel(name,
        read_compounds(book.sheet_by_name('metabolites')),
        read_reactions(book.sheet_by_name('reactions')))

    reaction_id, compound_name = model.check_reaction_compounds()
    if compound_name is not None:
        raise ParseError('Compound {}, {} not defined in compound table'.format(
            reaction_id, compound_name))

    return model


def read_modelseed(name, excel_path, ptt_file):
    """Read metabolic model for a ModelSEED model"""

    book = xlrd.open_workbook(excel_path)

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
    sheet = book.sheet_by_name('Genes')
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

    def read_compounds(sheet):
        for i in range(1, sheet.nrows):
            compound_id, name, alt_name, formula, charge, _ = (
                sheet.row_values(i, end_colx=6))

            if compound_id.strip() == '':
                continue

            charge = int(charge) if charge != '' else None

            entry = CompoundEntry(id=compound_id, name=name, formula=formula,
                                  formula_neutral=None, charge=charge,
                                  kegg=None, cas=None)
            yield compound_id, entry

    def read_reactions(sheet):
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

            entry = ReactionEntry(id=reaction_id, name=name, genes=genes,
                                  equation=equation, subsystem=None, ec=ec)
            yield reaction_id, entry

    model = MetabolicModel(
        name, read_compounds(book.sheet_by_name('Compounds')),
        read_reactions(book.sheet_by_name('Reactions')))

    reaction_id, compound_name = model.check_reaction_compounds()
    if compound_name is not None:
        raise ParseError('Compound {}, {} not defined in compound table'.format(
            reaction_id, compound_name))

    return model
