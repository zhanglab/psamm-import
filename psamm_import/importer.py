
"""Generic native model importer"""

import sys
import os
import errno
import argparse
import logging
from collections import OrderedDict

import yaml
from psamm.datasource import modelseed

from .datasource import Importer

logger = logging.getLogger(__name__)


# Define custom dict representers for YAML
# This allows reading/writing Python OrderedDicts in the correct order.
# See: https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())


def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))


def recursive_subclasses(cls):
    """Yield all subclasses of a class recursively"""
    for subclass in cls.__subclasses__():
        for subsubclass in recursive_subclasses(subclass):
            yield subsubclass
        yield subclass


def encode_utf8(s):
    if isinstance(s, unicode):
        return s.encode('utf-8')
    return s


def model_compounds(model):
    """Yield model compounds as YAML dicts"""

    for compound_id, compound in sorted(model.compounds.iteritems()):
        d = OrderedDict()
        d['id'] = encode_utf8(compound_id)

        name = compound.properties.get('name')
        if name is not None:
            d['name'] = encode_utf8(name)

        formula = compound.properties.get('formula')
        if formula is not None:
            d['formula'] = str(formula)

        formula_neutral = compound.properties.get('formula_neutral')
        if formula_neutral is not None:
            d['formula_neutral'] = str(formula_neutral)

        charge = compound.properties.get('charge')
        if charge is not None:
            d['charge'] = int(charge)

        kegg = compound.properties.get('kegg')
        if kegg is not None:
            d['kegg'] = encode_utf8(kegg)

        cas = compound.properties.get('cas')
        if cas is not None:
            d['cas'] = encode_utf8(cas)

        yield d


def model_reactions(model):
    """Yield model reactions as YAML dicts"""

    for reaction_id, reaction in sorted(model.reactions.iteritems()):
        d = OrderedDict()
        d['id'] = encode_utf8(reaction_id)

        if hasattr(reaction, 'name') and reaction.name is not None:
            d['name'] = encode_utf8(reaction.name)
        if hasattr(reaction, 'genes') and reaction.genes is not None:
            d['genes'] = [encode_utf8(g) for g in reaction.genes]
        if hasattr(reaction, 'equation') and reaction.equation is not None:
            d['equation'] = encode_utf8(modelseed.format_reaction(
                reaction.equation))
        if (hasattr(reaction, 'subsystem') and
                reaction.subsystem is not None):
            d['subsystem'] = encode_utf8(reaction.subsystem)
        if hasattr(reaction, 'ec') and reaction.ec is not None:
            d['ec'] = encode_utf8(reaction.ec)

        yield d


def main():
    parser = argparse.ArgumentParser(
        description='Import from external model formats')
    parser.add_argument('--source', metavar='path', default='.',
                        help='Source directory or file')
    parser.add_argument('--dest', metavar='path', default='.',
                        help='Destination directory (default is ".")')
    parser.add_argument('format', help='Format to import ("list" to see all)')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Discover all available model importers
    importers = {}
    for importer_class in recursive_subclasses(Importer):
        importer_name = getattr(importer_class, 'name', None)
        importer_title = getattr(importer_class, 'title', None)
        if (importer_name is not None and
                importer_title is not None and
                importer_name not in importers):
            importers[importer_name] = importer_class

    # Print list of importers
    if args.format in ('list', 'help'):
        print('Available importers:')
        if len(importers) == 0:
            logger.error('No importers found!')
        else:
            for name, importer_class in sorted(importers.iteritems(),
                                               key=lambda x: x[1].title):
                print('{:<10}  {}'.format(name, importer_class.title))
        sys.exit(0)

    if args.format not in importers:
        logger.error('Importer {} not found!'.format(args.format))
        logger.info('Use "list" to see available importers.')
        sys.exit(-1)

    importer = importers[args.format]()

    try:
        model = importer.import_model(args.source)
    except:
        logger.error('Failed to load model!', exc_info=True)
        importer.help()
        sys.exit(-1)

    model.print_summary()

    # Create destination directory if not exists
    dest = args.dest
    try:
        os.makedirs(dest)
    except OSError as e:
        if e.errno != errno.EEXIST or not os.path.isdir(dest):
            raise

    yaml.add_representer(OrderedDict, dict_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                         dict_constructor)

    yaml_args = {'default_flow_style': False,
                 'encoding': 'utf-8',
                 'allow_unicode': True,
                 'width': 79}

    with open(os.path.join(dest, 'compounds.yaml'), 'w+') as f:
        yaml.dump(list(model_compounds(model)), f, **yaml_args)

    with open(os.path.join(dest, 'reactions.yaml'), 'w+') as f:
        yaml.dump(list(model_reactions(model)), f, **yaml_args)

    model_d = OrderedDict([('name', model.name)])
    if model.biomass_reaction is not None:
        model_d['biomass'] = model.biomass_reaction
    model_d.update([
        ('compounds', [{'include': 'compounds.yaml'}]),
        ('reactions', [{'include': 'reactions.yaml'}])
    ])

    with open(os.path.join(dest, 'model.yaml'), 'w+') as f:
        yaml.dump(model_d, f, **yaml_args)
