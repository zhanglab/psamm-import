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

"""BiGG to native model importer."""

import sys
import os
import argparse
import logging
import json
from urllib import quote
import urllib2

from six import text_type

from .importer import write_yaml_model
from .datasource import cobrajson
from .util import mkdir_p
from .model import ParseError, ModelLoadError


logger = logging.getLogger(__name__)


def main():
    """Entry point for BiGG import program."""
    parser = argparse.ArgumentParser(
        description='Import from BiGG database')
    parser.add_argument('--dest', metavar='path', default='.',
                        help='Destination directory (default is ".")')
    parser.add_argument('--no-exchange', action='store_true',
                        help=('Disable importing exchange reactions as'
                              ' exchange compound file.'))
    parser.add_argument('--split-subsystem', action='store_true',
                        help='Enable splitting reaction files by subsystem')
    parser.add_argument('--force', action='store_true',
                        help='Enable overwriting model files')
    parser.add_argument('id', help='BiGG model to import ("list" to see all)')

    args = parser.parse_args()

    # Set up logging for the command line interface
    if 'PSAMM_DEBUG' in os.environ:
        level = getattr(logging, os.environ['PSAMM_DEBUG'].upper(), None)
        if level is not None:
            logging.basicConfig(level=level)
    else:
        logging.basicConfig(
            level=logging.INFO, format=u'%(levelname)s: %(message)s')

    # Print list of available models
    if args.id == 'list':
        print('Available models:')
        f = urllib2.urlopen('http://bigg.ucsd.edu/api/v2/models')
        results = json.load(f)['results']
        id_width = min(max(len(result['bigg_id']) for result in results), 16)
        for result in sorted(results, key=lambda x: x.get('organism')):
            print('{} {}'.format(
                result.get('bigg_id').ljust(id_width), result.get('organism')))
        sys.exit(0)

    importer = cobrajson.Importer()

    try:
        f = urllib2.urlopen(
            'http://bigg.ucsd.edu/api/v2/models/{}/download'.format(
                quote(args.id)))
        model = importer.import_model(f)
    except ModelLoadError as e:
        logger.error('Failed to load model!', exc_info=True)
        importer.help()
        parser.error(text_type(e))
    except ParseError as e:
        logger.error('Failed to parse model!', exc_info=True)
        logger.error(text_type(e))
        sys.exit(-1)

    model.print_summary()

    # Check if dest directory is empty. If we get an error assume that the
    # directory does not exist.
    dest_is_empty = False
    try:
        dest_is_empty = len(os.listdir(args.dest)) == 0
    except OSError:
        dest_is_empty = True

    if not dest_is_empty:
        if not args.force:
            logger.error('Destination directory is not empty. Use --force'
                         ' option to proceed anyway, overwriting any existing'
                         ' files in {}'.format(args.dest))
            return 1
        else:
            logger.warning('Destination directory is not empty, overwriting'
                           ' existing files in {}'.format(args.dest))

    # Create destination directory if not exists
    dest = args.dest
    mkdir_p(dest)

    convert_exchange = not args.no_exchange
    write_yaml_model(model, dest, convert_exchange=convert_exchange,
                     split_subsystem=args.split_subsystem)
