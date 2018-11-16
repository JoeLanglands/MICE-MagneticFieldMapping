import sys
import argparse

sys.path.insert(0, '..')
from micemag import core
from micemag.utils import default_grid

magDict = {'magnet': 'SSU', 'CC': {'I': 205.7, 'fb': None},\
           'M1': {'I': 0.0, 'fb': None}, 'M2': {'I': 0.0, 'fb': None}}

if __name__=='__main__':
    #core.buildG4BLfield(magDict, default_grid)
    core.showFBfield('ECENewData_Residual.dat', 'CC', 'SSU', nCores=5)
