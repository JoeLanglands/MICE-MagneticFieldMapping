import os
import pickle
import sys
import time

import numpy as np

import utils
from fbutils import applyfb as appFB
from fbutils import fbfit as fitFB
from fieldmanip.readData import readFile
from plotting import plots3d as p3d

#Simple functional interfaces to the fourier bessel fitting routines

def performFBfit(residField, coil, magnet, zmax=None, rmax=0.15, n=3, l=20, m=10,\
                 verbose=True, saveAs=None):
    if zmax==None:
        if coil=='ECE':
            zmax = 1.8
        else:
            zmax = 1.0

    if type(residField) == type('string'):
        fb_cls = fitFB.FBfitClass(readFile(os.path.join(utils.resid_field_path, residField)), \
                                  coil, magnet, zmax, rmax, n, l, m, verbose, saveAs)
    else:
        fb_cls = fitFB.FBfitClass(residField, coil, magnet, zmax, rmax, n, l, m, \
                                  verbose, saveAs)
    
    fb_cls.run()


def showFBfield(_residField, coil, magnet, fitDict=None, nCores=1):
    if type(_residField) == type('string'):
        residField = readFile(os.path.join(utils.resid_field_path, _residField))
    else:
        residField = _residField
    
    if fitDict == None:
        _fitDict = appFB.getDefaultFitDict(coil, magnet)
        
    fb_field = appFB.applyFB_field(residField, _fitDict, coil, magnet, FBonly=True, nCores=nCores)

    p3d.wireFrame(residField, fb_field)
    

def buildG4BLfield(magDict, gridDict, saveAs=None, FBonly=True):
    if saveAs == None:
        _date = time.localtime()
        saveAs = '%s_%s%s%s.table'%(magDict['magnet'], _date.tm_year, _date.tm_mon, _date.tm_mday)

    xNsteps = int((gridDict['x']['end'] + gridDict['x']['step'])/gridDict['x']['step'])
    xARR = np.linspace(gridDict['x']['start'], gridDict['x']['end'], xNsteps)
    
    yNsteps = int((gridDict['y']['end'] + gridDict['y']['step'])/gridDict['y']['step'])
    yARR = np.linspace(gridDict['y']['start'], gridDict['y']['end'], yNsteps)

    zNsteps = int((gridDict['z']['end'] + gridDict['z']['step'])/gridDict['z']['step'])
    zARR = np.linspace(gridDict['z']['start'], gridDict['z']['end'], zNsteps)

    scaleList = [' 1 X [1e3]\n', ' 2 Y [1e3]\n', ' 3 Z [1e3]\n', \
                 ' 4 BX [1e-3]\n', ' 5 BY [1e-3]\n', ' 6 BZ [1e-3]\n', ' 0\n']
    print 'Writing out %d field points'%(xNsteps*yNsteps*zNsteps)
    count = 1
    with open(os.path.join(utils.maus_field_path, saveAs), 'w') as _output:
        _output.write('\t%d\t%d\t%d\t1\n'%(xNsteps, yNsteps, zNsteps))
        for i in scaleList:
            _output.write(i)

        for _x in xARR:
            for _y in yARR:
                for _z in zARR:
                    if FBonly == True:
                        Bx, By, Bz = appFB.applyFB_grid(magDict, _x, _y, _z, 0, 0, 0)
                    elif FBonly == False:
                        #placeholder -- this needs the geofit interpolation added
                        Bx, By, Bz = appFB.applyFB_grid(magDict, x, y, z, 0, 0, 0) 
                    _output.write('{:.3f}\t{:.3f}\t{:.3f}\t{:.8f}\t{:.8f}\t{:.8f}\n'.format( \
                                    _x, _y,_z, Bx, By, Bz))
                    utils.progressBar(count, xNsteps*yNsteps*zNsteps)
                    count += 1
    print 'Finished! File can be found at %s'%os.path.join(utils.maus_field_path, saveAs)
                        
