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
# Copyright 2015  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""Importers for various Excel formats."""

import os
import re
import csv
import glob

import xlrd
from six import string_types

from psamm.datasource import native
from psamm.datasource.reaction import ReactionParser
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry,
                                    DictReactionEntry as ReactionEntry)
from psamm.datasource.context import FileMark, FilePathContext
from psamm.reaction import Reaction, Compound, Direction
from psamm.expression import boolean
from psamm.importer import Importer, ModelLoadError


class ImportiMA945(Importer):
    """Importer for iMA945 model."""

    name = 'iMA945'
    title = 'Salmonella enterica iMA945 (Excel format), AbuOun et al., 2009'

    filename = 'jbc.M109.005868-5.xls'

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.biomass_reaction = 'ST_biomass_core'
        model.extracellular_compartment = 'e'
        model.reactions.update(self._read_reactions())
        model.compounds.update(self._read_compounds())

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

            formula_neutral = self._try_parse_formula(
                compound_id, formula_neutral)

            # Skip formulas where the charge value is accidentally in the
            # formula column.
            if isinstance(formula, string_types):
                # Remove weird prefix. This could perhaps be charge values that
                # were accidentally put into the formula column.
                m = re.match(r'^(.*)-\d$', formula)
                if m:
                    formula = m.group(1)
                formula = self._try_parse_formula(compound_id, formula)
            else:
                formula = None

            kegg = None if kegg == '' else kegg

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, formula=formula,
                formula_neutral=formula_neutral,
                charge=charge, kegg=kegg), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<==>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows, parse_global=True)

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
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            genes = self._try_parse_gene_association(reaction_id, genes)

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name,
                genes=genes, equation=equation), filemark=filemark)


class ImportiRR1083(Importer):
    """Importer for iRR1083 model."""

    name = 'iRR1083'
    title = ('Salmonella enterica iRR1083 (Excel format),'
             ' Raghunathan et al., 2009')

    filename = '1752-0509-3-38-s1.xls'

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.extracellular_compartment = 'e'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

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

            formula_neutral = self._try_parse_formula(
                compound_id, formula_neutral)

            kegg = None if kegg == '' else kegg

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name,
                formula=formula_neutral,
                formula_neutral=formula_neutral,
                charge=charge, kegg=kegg), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<==>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows, parse_global=True)

        sheet = self._book.sheet_by_name('Gene Protein Reaction iRR1083')
        for i in range(3, sheet.nrows):
            genes, protein, reaction_id, name, equation, subsystem = (
                sheet.row_values(i, end_colx=6))

            if reaction_id.strip() == '':
                continue

            genes = self._try_parse_gene_association(reaction_id, genes)
            protein = None if protein == '' else protein
            name = None if name == '' else name

            if equation != '':
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            subsystem = None if subsystem == '' else subsystem

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation), filemark=filemark)


class ImportiJO1366(Importer):
    """Importer for iJO1366 model."""

    name = 'iJO1366'
    title = ('Escerichia coli iJO1366 (Excel format),'
             ' Orth et al., 2011')

    filename = 'inline-supplementary-material-2.xls'

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M'
        model.extracellular_compartment = 'e'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

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

            formula_neutral = self._try_parse_formula(
                compound_id, formula_neutral)
            formula = self._try_parse_formula(compound_id, formula)

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            kegg = None if kegg.strip() == '' else kegg
            cas = None if cas.strip() == '' else cas

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name,
                formula=formula,
                formula_neutral=formula_neutral,
                charge=charge, kegg=kegg, cas=cas), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('->', Direction.Forward),
            ('<=>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows)

        sheet = self._book.sheet_by_name('Table 2')
        for i in range(1, sheet.nrows):
            (reaction_id, name, equation, _, genes, _, subsystem, ec,
                reversible) = sheet.row_values(i, end_colx=9)

            if reaction_id.strip() == '':
                continue

            genes = self._try_parse_gene_association(reaction_id, genes)
            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            ec = None if ec.strip() == '' else ec

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, subsystem=subsystem,
                ec=ec), filemark=filemark)


class EColiTextbookImport(Importer):
    """Importer for E. coli core textbook model."""

    name = 'EColi_textbook'
    title = ('Escerichia coli Textbook (core) model (Excel format),'
             ' Orth et al., 2010')

    filename = 'ecoli_core_model.xls'

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.extracellular_compartment = 'e'
        model.reactions.update(self._read_reactions())
        model.compounds.update(self._read_compounds())

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

            formula_neutral = self._try_parse_formula(
                compound_id, formula_neutral)
            formula = self._try_parse_formula(compound_id, formula)

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            cas = None if cas.strip() == '' or cas == 'None' else cas
            kegg = None if kegg.strip() == '' else kegg

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name,
                formula=formula,
                formula_neutral=formula_neutral,
                charge=charge, kegg=kegg, cas=cas), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<==>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows, parse_global=True)

        sheet = self._book.sheet_by_name('reactions')
        for i in range(1, sheet.nrows):
            (reaction_id, name, equation, subsystem, ec, _, _, _, _, _,
                genes) = sheet.row_values(i, end_colx=11)

            if reaction_id.strip() == '':
                continue

            genes = self._try_parse_gene_association(reaction_id, genes)
            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            ec = None if ec.strip() == '' else ec

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes, equation=equation,
                subsystem=subsystem, ec=ec), filemark=filemark)


class ImportSTMv1_0(Importer):  # noqa
    """Importer for STM_v1.0 model."""

    name = 'STM_v1.0'
    title = ('Salmonella enterica STM_v1.0 (Excel format),'
             ' Thiele et al., 2011')

    filename = '1752-0509-5-8-s1.xlsx'

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.biomass_reaction = 'biomass_iRR1083_metals'
        model.extracellular_compartment = 'e'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('SI Tables - S2b - Metabolites')
        for i in range(2, sheet.nrows):
            _, compound_id, name, formula, charge, _, kegg, pubchem, chebi = (
                sheet.row_values(i, end_colx=9))

            if compound_id.strip() == '':
                continue

            name = None if name.strip() == '' else name.strip()
            formula = self._try_parse_formula(compound_id, formula)

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            kegg = None if kegg.strip() == '' else kegg

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name,
                formula=formula, charge=charge, kegg=kegg), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<=>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows)

        sheet = self._book.sheet_by_name('SI Tables - S2a - Reactions')
        for i in range(4, sheet.nrows):
            reaction_id, name, equation, genes, _, subsystem = (
                sheet.row_values(i, end_colx=6))

            if reaction_id.strip() == '':
                continue

            name = None if name.strip() == '' else name.strip()

            if equation.strip() != '':
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            genes = self._try_parse_gene_association(reaction_id, genes)

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name,
                genes=genes, equation=equation,
                subsystem=subsystem), filemark=filemark)


class ImportiJN746(Importer):
    """Importer for iJN746 model."""

    name = 'iJN746'
    title = ('Pseudomonas putida iJN746 (Excel format),'
             ' Nogales et al., 2011')

    filenames = ('1752-0509-2-79-s8.xls',
                 '1752-0509-2-79-s9.xls')

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:')
        for filename in self.filenames:
            print('- {}'.format(filename))

    def import_model(self, source):
        """Import and return model instance."""
        if not os.path.isdir(source):
            raise ModelLoadError('Source must be a directory')

        self._compound_context = FilePathContext(
            os.path.join(source, self.filenames[0]))
        self._reaction_context = FilePathContext(
            os.path.join(source, self.filenames[1]))

        self._compound_book = xlrd.open_workbook(
            self._compound_context.filepath)
        self._reaction_book = xlrd.open_workbook(
            self._reaction_context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.extracellular_compartment = 'e'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._compound_book.sheet_by_name('Additional file 8')
        for i in range(1, sheet.nrows):
            (compound_id, name, formula, charge, cas, formula_neutral, _,
                kegg) = sheet.row_values(i, end_colx=8)

            if compound_id.strip() == '':
                continue

            name = None if name.strip() == '' else name

            formula_neutral = self._try_parse_formula(
                compound_id, formula_neutral)
            formula = self._try_parse_formula(compound_id, formula)

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            kegg = None if kegg.strip() == '' else kegg

            if isinstance(cas, basestring) and cas.strip() in ('', 'None'):
                cas = None
            else:
                cas = str(cas)

            filemark = FileMark(self._compound_context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, formula=formula,
                formula_neutral=formula_neutral,
                charge=charge, kegg=kegg, cas=cas), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<==>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows, parse_global=True)

        sheet = self._reaction_book.sheet_by_name('Additional file 9')
        for i in range(1, sheet.nrows):
            reaction_id, name, equation, subsystem, ec, _, genes = (
                sheet.row_values(i, end_colx=7))

            if reaction_id.strip() == '':
                continue

            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            ec = None if ec.strip() == '' else ec
            genes = self._try_parse_gene_association(reaction_id, genes)

            filemark = FileMark(self._reaction_context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name,
                genes=genes, equation=equation,
                subsystem=subsystem, ec=ec), filemark=filemark)


class ImportiJP815(Importer):
    """Importer for iJP815 model."""

    name = 'iJP815'
    title = ('Pseudomonas putida iJP815 (Excel format),'
             ' Puchalka et al., 2008')

    filename = 'journal.pcbi.1000210.s011.XLS'

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

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

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, kegg=kegg), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<==>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows)

        sheet = self._book.sheet_by_name('Reactions')
        for i in range(1, sheet.nrows):
            (reaction_id, name, equation, _, _, _, _, _, _, subsystem,
                genes) = sheet.row_values(i, end_colx=11)

            if reaction_id.strip() == '':
                continue

            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)

                # Rebuild reaction with compartment information
                def translate(c, v):
                    compartment = 'e' if c.name[0] == 'E' else None
                    return Compound(c.name[1:], compartment=compartment), v

                equation = Reaction(
                    equation.direction,
                    (translate(c, v) for c, v in equation.compounds))
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            genes = self._try_parse_gene_association(reaction_id, genes)

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, subsystem=subsystem), filemark=filemark)


class ImportiSyn731(Importer):
    """Importer for iSyn731."""

    name = 'iSyn731'
    title = ('Synechocystis sp. PCC 6803 iSyn731 (Excel format),'
             ' Saha et al., 2012')

    filename = 'journal.pone.0048285.s001.XLSX'

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.biomass_reaction = 'Biomass_Hetero'
        model.extracellular_compartment = 'e'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('Metabolites')
        for i in range(1, sheet.nrows):
            compound_id, name, formula, charge, kegg = (
                sheet.row_values(i))

            if compound_id.strip() == '':
                continue

            name = None if name == '' else name
            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            formula = formula.strip()
            if formula not in ('', '-', 'noformula'):
                formula = self._try_parse_formula(compound_id, formula)
            else:
                formula = None

            if kegg != 0 and kegg.strip() != '':
                kegg = kegg.split('|')
                kegg = kegg if len(kegg) > 1 else kegg[0]
            else:
                kegg = None

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, formula=formula,
                charge=charge, kegg=kegg), filemark=filemark)

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
                # check that this works correctly. should substitute the => for
                # a space then =>. the double spaces should be ignored though.
                equation = re.sub(r'\s*\+\s*', ' + ', equation)
                equation = re.sub(r'\|\[(\w)\]', r'[\1]|', equation)
                equation = equation.replace('||', '|')
                equation = self._try_parse_reaction(reaction_id, equation)
            else:
                equation = None

            genes = self._try_parse_gene_association(reaction_id, genes)
            ec = ec if ec.strip() != '' and ec != 'Undetermined' else None

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, subsystem=subsystem,
                ec=ec), filemark=filemark)


class ImportiCce806(Importer):
    """Importer for iCce806 model."""

    name = 'iCce806'
    title = ('Cyanothece sp. ATCC 51142 iCce806 (Excel format),'
             ' Vu et al., 2012')

    filenames = ('journal.pcbi.1002460.s005.XLSX',
                 'journal.pcbi.1002460.s006.XLSX')

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:')
        for filename in self.filenames:
            print('- {}'.format(filename))

    def import_model(self, source):
        """Import and return model instance."""
        if not os.path.isdir(source):
            raise ModelLoadError('Source must be a directory')

        self._compound_context = FilePathContext(
            os.path.join(source, self.filenames[1]))
        self._reaction_context = FilePathContext(
            os.path.join(source, self.filenames[0]))

        self._compound_book = xlrd.open_workbook(
            self._compound_context.filepath)
        self._reaction_book = xlrd.open_workbook(
            self._reaction_context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.biomass_reaction = 'CyanoBM (average)'
        model.extracellular_compartment = 'e'
        model.reactions.update(self._read_reactions())
        model.compounds.update(self._read_compounds())

        return model

    def _read_compounds(self):
        sheet = self._compound_book.sheet_by_name('Table S2')
        for i in range(2, sheet.nrows):
            (compound_id, name, formula, charge, cas, formula_neutral, _,
                kegg) = sheet.row_values(i)

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

            m = re.match(r'^(.*)-\d$', formula)
            if m:
                formula = m.group(1)
            formula = self._try_parse_formula(compound_id, formula)

            kegg = None if kegg == '' else kegg

            if cas.strip() != '' and cas.strip() != 'None':
                cas = cas.strip()
            else:
                cas = None

            filemark = FileMark(self._compound_context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, formula=formula,
                formula_neutral=formula, charge=charge,
                kegg=kegg, cas=cas), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<==>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows, parse_global=True)

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
                equation = re.sub(r'\s*\+\s*', ' + ', equation)
                equation = re.sub(r'\s*-->\s*', ' --> ', equation)
                equation = re.sub(r'\s*<==>\s*', ' <==> ', equation)
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            genes = self._try_parse_gene_association(reaction_id, genes)

            m = re.match(r'EC-(.*)', ec)
            if m:
                ec = m.group(1)
                if ec == 'Undetermined':
                    ec = None
            else:
                ec = None

            filemark = FileMark(self._reaction_context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, subsystem=subsystem,
                ec=ec), filemark=filemark)


class ImportGSMN_TB(Importer):  # noqa
    """Importer for GSMN-TB model."""

    name = 'GSMN-TB'
    title = ('Mycobacterium tuberculosis GSMN-TB (Excel format),'
             ' Beste et al., 2007')

    filenames = ('gb-2007-8-5-r89-s4.xls',
                 'gb-2007-8-5-r89-s6.xls')

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:')
        for filename in self.filenames:
            print('- {}'.format(filename))

    def import_model(self, source):
        """Import and return model instance."""
        if not os.path.isdir(source):
            raise ModelLoadError('Source must be a directory')

        self._compound_context = FilePathContext(
            os.path.join(source, self.filenames[1]))
        self._reaction_context = FilePathContext(
            os.path.join(source, self.filenames[0]))

        self._compound_book = xlrd.open_workbook(
            self._compound_context.filepath)
        self._reaction_book = xlrd.open_workbook(
            self._reaction_context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

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

            filemark = FileMark(self._compound_context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name), filemark=filemark)

        def create_missing(compound_id, name=None):
            if name is None:
                name = compound_id
            return CompoundEntry(dict(id=compound_id, name=name))

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
        arrows = (
            ('->', Direction.Forward),
            ('=', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows)

        sheet = self._reaction_book.sheet_by_name('File 4')
        for i in range(4, sheet.nrows):
            (reaction_id, equation, fluxbound, _, ec, genes, name,
                subsystem) = sheet.row_values(i, end_colx=8)

            if reaction_id.startswith('%') or reaction_id.strip() == '':
                continue
            genes = self._try_parse_gene_association(reaction_id, genes)

            name = None if name.strip() == '' else name
            if equation.strip() != '':
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
                rdir = Direction.Both if fluxbound != 0 else Direction.Forward
                equation = Reaction(rdir, equation.compounds)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem
            ec = None if ec.strip() == '' else ec

            filemark = FileMark(self._reaction_context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, subsystem=subsystem,
                ec=ec), filemark=filemark)


class ImportiNJ661(Importer):
    """Importer for iNJ661 model."""

    name = 'iNJ661'
    title = ('Mycobacterium tuberculosis iNJ661 (Excel format),'
             ' Jamshidi et al., 2007')

    filename = '1752-0509-1-26-s5.xls'

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model(self, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = self.title
        model.extracellular_compartment = 'e'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('metabolites')
        for i in range(1, sheet.nrows):
            compound_id, name, formula, charge = sheet.row_values(
                i, end_colx=4)

            if compound_id.strip() == '':
                continue

            name = name if name.strip() != '' else None
            formula = self._try_parse_formula(compound_id, formula)

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, formula=formula,
                charge=charge), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<==>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows, parse_global=True)

        sheet = self._book.sheet_by_name('iNJ661')
        for i in range(5, sheet.nrows):
            reaction_id, name, equation, _, subsystem, _, genes = (
                sheet.row_values(i, end_colx=7))

            if reaction_id.strip() == '':
                continue

            # TODO model uses an alternative gene association format
            if genes.strip() != '':
                genes = frozenset(m.group(0)
                                  for m in re.finditer(r'Rv\w+', genes))
            else:
                genes = None

            name = None if name.strip() == '' else name

            if equation.strip() != '':
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, subsystem=subsystem), filemark=filemark)


class ImportGenericiNJ661mv(Importer):
    """Importer for the models iNJ661m and iNJ661v.

    For models of Mycobacterium tuberculosis iNJ661m/v (Excel format),
    Fang et al., 2010.
    """

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model_named(self, name, source):
        """Import and return model instance with the given name."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)

        model = native.NativeModel()
        model.name = name
        model.biomass_reaction = 'biomass_Mtb_9_60atp_test_NOF'
        model.extracellular_compartment = 'e'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('metabolites')
        for i in range(1, sheet.nrows):
            compound_id, name, formula = sheet.row_values(i, end_colx=3)

            # Skip compartmentalized compounds
            compound_id = re.match(r'^(.*)\[.\]$', compound_id).group(1)
            name = name if name.strip() != '' else None

            formula = self._try_parse_formula(compound_id, formula)

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, formula=formula), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('->', Direction.Forward),
            ('<=>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows)

        sheet = self._book.sheet_by_name('reactions')
        for i in range(1, sheet.nrows):
            reaction_id, name, equation, genes, _, subsystem = (
                sheet.row_values(i, end_colx=6))

            if reaction_id.strip() == '':
                continue

            genes = self._try_parse_gene_association(reaction_id, genes)
            name = None if name.strip() == '' else name.strip()

            # Biomass reaction is not specified in this table
            if (equation.strip() != '' and
                    reaction_id != 'biomass_Mtb_9_60atp_test_NOF'):
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
            else:
                equation = None

            subsystem = None if subsystem.strip() == '' else subsystem

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, subsystem=subsystem), filemark=filemark)


class ImportiNJ661m(ImportGenericiNJ661mv):
    """Importer for iNJ661m model."""

    name = 'inj661m'
    title = ('Mycobacterium tuberculosis iNJ661m (Excel format),'
             ' Fang et al., 2010')
    filename = '1752-0509-4-160-s3.xls'

    def import_model(self, source):
        """Import and return model instance."""
        return self.import_model_named(self.title, source)


class ImportiNJ661v(ImportGenericiNJ661mv):
    """Importer for iNJ661v model."""

    name = 'inj661v'
    title = ('Mycobacterium tuberculosis iNJ661v (Excel format),'
             ' Fang et al., 2010')
    filename = '1752-0509-4-160-s5.xls'

    def import_model(self, source):
        """Import and return model instance."""
        return self.import_model_named(self.title, source)


class ImportShewanellaOng(Importer):
    """Generic importer for four models published in Ong et al., 2014.

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
        """Print import help text."""
        print('Source must contain the model definition in Excel format.\n'
              'Expected files in source directory:\n'
              '- {}'.format(self.filename))

    def import_model_named(self, name, col_index, source):
        """Import and return model instance."""
        context = FilePathContext(source)
        if os.path.isdir(context.filepath):
            context = FilePathContext(os.path.join(source, self.filename))

        self._context = context
        self._book = xlrd.open_workbook(context.filepath)
        self._col_index = col_index

        model = native.NativeModel()
        model.name = name
        model.biomass_reaction = self.biomass_names[col_index]
        model.extracellular_compartment = 'e'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions())

        return model

    def _read_compounds(self):
        sheet = self._book.sheet_by_name('S3-Metabolites')
        for i in range(1, sheet.nrows):
            (compound_id, _, _, _, _, _, name, formula_neutral, formula,
                charge, _, kegg, cas) = sheet.row_values(i, end_colx=13)

            # Remove compartmentalization of compounds
            compound_id = re.match(
                r'^(.*)\[.\]$', compound_id).group(1).lower()
            name = name if name.strip() != '' else None

            formula_neutral = self._try_parse_formula(
                compound_id, formula_neutral)
            formula = self._try_parse_formula(compound_id, formula)

            try:
                charge = None if charge == '' else int(charge)
            except ValueError:
                charge = None

            kegg = None if kegg.strip() == '' else kegg

            if isinstance(cas, basestring) and cas.strip() in ('', 'None'):
                cas = None
            else:
                cas = str(cas)

            filemark = FileMark(self._context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, formula=formula,
                formula_neutral=formula_neutral, charge=charge,
                kegg=kegg, cas=cas), filemark=filemark)

    def _read_reactions(self):
        arrows = (
            ('-->', Direction.Forward),
            ('<==>', Direction.Both)
        )
        parser = ReactionParser(arrows=arrows, parse_global=True)

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
                i, start_colx=8, end_colx=12)[self._col_index].strip()

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
                equation = re.sub(r'\s*\+\s*', ' + ', equation)
                equation = self._try_parse_reaction(
                    reaction_id, equation, parser=parser.parse)
                equation = equation.translated_compounds(translate)
            else:
                equation = None

            genes = self._try_parse_gene_association(reaction_id, model_genes)
            subsystem = sheet.cell_value(i, 18)

            name = None if name.strip() == '' else name.strip()
            subsystem = None if subsystem.strip() == '' else subsystem

            filemark = FileMark(self._context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, subsystem=subsystem), filemark=filemark)


class ImportiMR1_799(ImportShewanellaOng):  # noqa
    """Importer for iMR_799 model."""

    name = 'imr1_799'
    title = ('Shewanella oneidensis MR-1 iMR1_799 (Excel format),'
             ' Ong et al., 2014')

    def import_model(self, source):
        """Import and return model instance."""
        return self.import_model_named(self.title, 0, source)


class ImportiMR4_812(ImportShewanellaOng):  # noqa
    """Importer for iMR4_812 model."""

    name = 'imr4_812'
    title = ('Shewanella sp. MR-4 iMR4_812 (Excel format),'
             ' Ong et al., 2014')

    def import_model(self, source):
        """Import and return model instance."""
        return self.import_model_named(self.title, 1, source)


class ImportiW3181_789(ImportShewanellaOng):  # noqa
    """Importer for iW3181_789 model."""

    name = 'iw3181_789'
    title = ('Shewanella sp. W3-18-1 iW3181_789 (Excel format),'
             ' Ong et al., 2014')

    def import_model(self, source):
        """Import and return model instance."""
        return self.import_model_named(self.title, 2, source)


class ImportiOS217_672(ImportShewanellaOng):  # noqa
    """Importer for iOS217_672 model."""

    name = 'ios217_672'
    title = ('Shewanella denitrificans OS217 iOS217_672 (Excel format),'
             ' Ong et al., 2014')

    def import_model(self, source):
        """Import and return model instance."""
        return self.import_model_named(self.title, 3, source)


class ImportModelSEED(Importer):
    """Read metabolic model for a ModelSEED model."""

    name = 'ModelSEED'
    title = 'ModelSEED model (Excel format)'
    generic = True

    def help(self):
        """Print importer help text."""
        print('Source must contain the model definition in Excel format\n'
              ' and a PTT file for mapping PEG gene names.'
              'Expected files in source directory:\n'
              '- Seed*.xls\n'
              '- NC_*.ptt')

    def import_model(self, source):
        """Import and return model instance."""
        if not os.path.isdir(source):
            raise ModelLoadError('Source must be a directory')

        excel_sources = glob.glob(os.path.join(source, 'Seed*.xls'))
        if len(excel_sources) == 0:
            raise ModelLoadError('No .xls file found in source directory')
        elif len(excel_sources) > 1:
            raise ModelLoadError(
                'More than one .xls file found in source directory')

        self._excel_context = FilePathContext(excel_sources[0])

        ptt_sources = glob.glob(os.path.join(source, '*.ptt'))
        if len(ptt_sources) == 0:
            raise ModelLoadError('No .ptt file found in source directory')
        elif len(ptt_sources) > 1:
            raise ModelLoadError(
                'More than one .ptt file found in source directory')

        self._book = xlrd.open_workbook(self._excel_context.filepath)

        with open(ptt_sources[0], 'r') as ptt_file:
            # Read mapping from location to gene ID from PTT file
            location_mapping = {}
            for i in range(3):
                ptt_file.readline()  # Skip headers
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

        model = native.NativeModel()
        model.name = 'ModelSEED model'
        model.compounds.update(self._read_compounds())
        model.reactions.update(self._read_reactions(peg_mapping))

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

            filemark = FileMark(self._excel_context, i, None)
            yield CompoundEntry(dict(
                id=compound_id, name=name, formula=formula,
                charge=charge), filemark=filemark)

    def _read_reactions(self, peg_mapping):
        sheet = self._book.sheet_by_name('Reactions')
        for i in range(1, sheet.nrows):
            reaction_id, name, equation, _, ec_list, _, _, pegs = (
                sheet.row_values(i, end_colx=8))

            if reaction_id.strip() == '':
                continue

            name = name if name.strip() != '' else None

            if equation != '' and 'NONE' not in equation:
                equation = self._try_parse_reaction(reaction_id, equation)
            else:
                continue

            def translate_pegs(variable):
                if variable.symbol in peg_mapping:
                    return boolean.Variable(peg_mapping[variable.symbol])
                return True

            pegs = self._try_parse_gene_association(reaction_id, pegs)
            genes = None
            if isinstance(pegs, boolean.Expression):
                genes = pegs.substitute(translate_pegs)
                if genes.has_value():
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

            filemark = FileMark(self._excel_context, i, None)
            yield ReactionEntry(dict(
                id=reaction_id, name=name, genes=genes,
                equation=equation, ec=ec), filemark=filemark)
