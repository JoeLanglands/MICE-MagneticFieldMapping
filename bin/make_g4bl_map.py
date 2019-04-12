#!/usr/bin/env python2

"""Main script that generates g4blgrid fields from command line inputs.

   Can be run as executable after using $chmod +x make_g4bl_map.py.  Takes command line arguments.
   Example:
       $./make_g4bl_grid.py SSU -M1 185.0 -CC 205.7
   This would output the field from SSU with M1 set at 185.0 A, E1-CC-E2 set at 205.7 A and M2 at 0 A.
   The -h and --help options display the other options.

   Todo:
       *Needs argparse + dictionary support for using different fit dictionaries for the FB terms.
        Although this will be seldom used.
       *Needs support for the geometrical fit stuff when the code exists for that
       *Needs support for specifying the g4bl grid output since it only uses the default right now
"""

import sys
import argparse

sys.path.insert(0, '..') #Because setup.py does not work right now
from micemag import core
import micemag.utils as utils


parser = argparse.ArgumentParser(description='Program to make magnetic field in g4bl format for MAUS.')
parser.add_argument('magnet', help='Desired magnet: ssu/SSU or ssd/SSD')
parser.add_argument('-M1', '--Match1', help='Desired current of M1 in Amperes. Defaults to 0', \
                    type=float, default=0.0)

parser.add_argument('-M2', '--Match2', help='Desired current of M2 in Amperes. Defaults to 0', \
                    type=float, default=0.0)

parser.add_argument('-CC', '--Centre', help='Desired current of E1-CC-E2 in Amperes. Defaults to 0', \
                    type=float, default=0.0)

parser.add_argument('-FB', '--FBonly', help='Use this flag if you only want the FB contributions\n', \
                    action='store_true')

parser.add_argument('-t', '--time', help='Show execution time of program', action='store_true')

parser.add_argument('-s', '--saveas', help='Name to save output file as', default=None)

def buildDicts(args):
    if args.magnet not in ['ssu', 'SSU', 'ssd', 'SSD']:
        print 'Unrecognised magnet parameter. Please only use: ssu, SSU, ssd, SSD.'
        sys.exit()
    if args.Match1 == 0.0 and args.Match2 == 0.0 and args.Centre == 0.0:
        print 'You have not specified any currents!'
        sys.exit()  

    magDict = {'magnet': args.magnet.upper(), \
               'CC': {'I': args.Centre, 'fb': None}, \
               'M1': {'I': args.Match1, 'fb': None}, \
               'M2': {'I': args.Match2, 'fb': None}}

    return magDict
    
if __name__=='__main__':
    args = parser.parse_args()

    if args.time:
        @utils.timeit
        def make_g4bl_map_timed(args):
            core.buildG4BLfield(buildDicts(args),utils.default_grid, saveAs=args.saveas, \
                                FBonly=args.FBonly)
            return
        make_g4bl_map_timed(args)
    else:
        core.buildG4BLfield(buildDicts(args), utils.default_grid, saveAs=args.saveas, \
                            FBonly=args.FBonly)
    
    
