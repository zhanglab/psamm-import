#!/usr/bin/env python

import sys
import os
import errno
import argparse
import logging
from collections import OrderedDict

import yaml
from metnet.datasource import modelseed

from datasource import Importer

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


if __name__ == '__main__':
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
        logger.error('Importer {} not found!')
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

    def encode_utf8(s):
        if isinstance(s, unicode):
            return s.encode('utf-8')
        return s

    def model_compounds():
        for compound_id, compound in sorted(model.compounds.iteritems()):
            d = OrderedDict((key, value)
                            for key, value in compound._asdict().iteritems()
                            if value is not None)
            #d = {key: value for key, value in compound._asdict().iteritems()
            #     if value is not None}

            d['id'] = encode_utf8(d['id'])

            if 'name' in d:
                d['name'] = encode_utf8(d['name'])
            if 'formula' in d:
                d['formula'] = str(d['formula'])
            if 'formula_neutral' in d:
                d['formula_neutral'] = str(d['formula_neutral'])
            if 'kegg' in d:
                d['kegg'] = encode_utf8(d['kegg'])
            if 'cas' in d:
                d['cas'] = encode_utf8(d['cas'])

            yield d

    def model_reactions():
        for reaction_id, reaction in sorted(model.reactions.iteritems()):
            d = OrderedDict((key, value)
                            for key, value in reaction._asdict().iteritems()
                            if value is not None)

            d['id'] = encode_utf8(d['id'])

            if 'name' in d:
                d['name'] = encode_utf8(d['name'])
            if 'genes' in d:
                d['genes'] = [encode_utf8(g) for g in d['genes']]
            if 'equation' in d:
                d['equation'] = encode_utf8(modelseed.format_reaction(
                    d['equation']))
            if 'subsystem' in d:
                d['subsystem'] = encode_utf8(d['subsystem'])
            if 'ec' in d:
                d['ec'] = encode_utf8(d['ec'])

            yield d

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
                 'allow_unicode': True}

    with open(os.path.join(dest, 'compounds.yaml'), 'w+') as f:
        yaml.dump(list(model_compounds()), f, **yaml_args)

    with open(os.path.join(dest, 'reactions.yaml'), 'w+') as f:
        yaml.dump(list(model_reactions()), f, **yaml_args)

    with open(os.path.join(dest, 'model.yaml'), 'w+') as f:
        yaml.dump(OrderedDict((
            ('name', model.name),
            ('compounds', [{'include': 'compounds.yaml'}]),
            ('reactions', [{'include': 'reactions.yaml'}])
        )), f, **yaml_args)
