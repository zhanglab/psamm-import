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

"""BiGG to native model importer."""

import sys
import argparse
import logging
import json
from urllib import quote
import urllib2

from .importer import write_yaml_model
from .datasource import cobrajson
from .util import mkdir_p


logger = logging.getLogger(__name__)


def main():
    """Entry point for BiGG import program."""
    parser = argparse.ArgumentParser(
        description='Import from BiGG database')
    parser.add_argument('--dest', metavar='path', default='.',
                        help='Destination directory (default is ".")')
    parser.add_argument('--no-medium', action='store_true',
                        help='Disable importing exchange reactions as medium')
    parser.add_argument('id', help='BiGG model to import ("list" to see all)')

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Print list of importers
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
    except:
        logger.error('Failed to load model!', exc_info=True)
        sys.exit(-1)

    model.print_summary()

    # Create destination directory if not exists
    dest = args.dest
    mkdir_p(dest)

    write_yaml_model(model, dest, convert_medium=not args.no_medium)
