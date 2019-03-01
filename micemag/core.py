import os
import pickle
import sys
import time

import numpy as np

import matplotlib.pyplot as plt

import utils
from fbutils import applyfb as appFB
from fbutils import fbfit as fitFB
from fieldmanip.readData import readFile
from fieldmanip import fieldManipulation as fm
from fieldmanip import polarMeasurement as rphiz
from plotting import plots3d as p3d
from makefields import mkfieldclass as mkfield
from geofit import geofit
from geofit import coilfit

"""
This core module as the name suggests contains a few functions that could be considered as core
features of this package.  Everything that you definitely would want to do is defined as a function
here.

"""

def performFBfit(residField, magnet, coil, zmax=None, rmax=0.15, n=3, l=20, m=10,\
                 verbose=True, saveAs=None):
    if zmax==None:
        if coil in ['CC', 'ECE']:
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


def showFBfield(_residField, magnet, coil, fitDict=None, nCores=1):
    if type(_residField) == type('string'):
        residField = readFile(os.path.join(utils.resid_field_path, _residField))
    else:
        residField = _residField
    
    if fitDict == None:
        _fitDict = appFB.getDefaultFitDict(coil, magnet)
    else:
        with (os.path.join(utils.fb_pickle_path, fitDict), 'rb') as _pickle:
            _fitDict = pickle.load(_pickle)
        
    fb_field = appFB.applyFB_field(residField, _fitDict, coil, magnet, FBonly=True, nCores=nCores)

    p3d.wireFrame(residField, fb_field)
    

def buildG4BLfield(magDict, gridDict, saveAs=None, FBonly=False, coil=True):
    """Builds a magnetic field of SSU/SSD and prints it out to a .table file in g4blgrid format.

    Args:
        magDict  (dict): Dictionary containing magnet, coil currents and custom fitDict paths.
                         If fitDict paths are not specified it pulls the default ones.
        gridDict (dict): Dictionary containing information about the grid in which to calculate
                         the field over. 
        saveAs   (str):  Name that the user wishes to call the outputted field (no need to supply
                         full path). If None (default value), the magnet name + todays date is used.
        FBonly (bool):   When True: calculate only FB terms.  When False: calculate geofit+FB terms,
                         i.e the full model field is output.
        coil (bool):     When true, the full field is calculated from the coil fit model. If false,
                         the geometrical fit model is used instead.

    Returns:
        Doesn't return anything.  The outputted field is saved at data/MAUS/saveAs.table.

    Todo:
        *The scaleList part could change? May need support so that it can be adjusted by the user    
    
    """
    print 'Calculating field map for magnet:', magDict['magnet']
    print 'With currents:'
    print '\n\t M1  -> %.2f A\n\t M2  -> %.2f A\n\t ECE -> %.2f A\n'%(magDict['M1']['I'], \
                                                                      magDict['M2']['I'], \
                                                                      magDict['CC']['I'])
    if FBonly == False and coil == True:
        coilfit_calc = get_coilfit_class(magDict)
        
    print 'This could take a while...'
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

    
    start_time = time.time()
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
                        _Bx, _By, _Bz = coilfit_calc.calc_full_field_at_point_xyz(_x, _y, _z)
                        Bx, By, Bz = appFB.applyFB_grid(magDict, _x, _y, _z, _Bx, _By, _Bz)
                    _output.write('{:.3f}\t{:.3f}\t{:.3f}\t{:.8f}\t{:.8f}\t{:.8f}\n'.format( \
                                    _x, _y,_z, Bx, By, Bz))
                    utils.progressBar(count, xNsteps*yNsteps*zNsteps, start_time, time.time())
                    count += 1
                        
    print 'Finished! File can be found at %s'%os.path.join(utils.maus_field_path, saveAs)

    
    
def perform_coil_fit(magnet, coil, FBfit=False, makeresid=True, save_as=None, verbose=True):
    if magnet.upper() not in ['SSU', 'SSD']:
        print 'Magnet unrecognised - please use SSU or SSD'
        return
    if coil.upper() not in ['M1', 'M2', 'CC', 'ECE']:
        print 'Coil unrecognised - please use M1, M2, CC or ECE'
        print '\tN.B You can not fit to the end coils individually, only to E1-CC-E2'
        return
    if coil.upper() == 'CC':
        coil = 'ECE'

    if save_as == None:
        save_str = os.path.join(utils.geofit_field_path, magnet.upper() + '_' + coil.upper() \
                                + '_coilfit_default.pickle')
    else:
        save_str = os.path.join(utils.geofit_field_path, save_as)

    if coil.upper() in ['M1', 'M2']:
        print 'Performing coil fit on', magnet.upper(), coil.upper()
        if utils.coil_datacards[magnet.upper()][coil.upper()]['30A_data'] == None:
            print 'No data to fit to for this magnet!'
            return
        
        _centre = utils.centres_dict[magnet.upper()]['mapper'][coil.upper()]

        _field = readFile(utils.coil_datacards[magnet.upper()][coil.upper()]['30A_data'])

        if magnet.upper() == 'SSD':
            _field = fm.flip_SSD_data(_field)
        
        coilFitClass = coilfit.CoilFitClass(utils.coil_datacards[magnet.upper()][coil.upper()], \
                                              _field, _centre)
        fitDict = coilFitClass.run()
        print 'Finished with parameters: '
        for key, value in fitDict.iteritems():
            print key, value
        print 'Saved fit parameters at: ', save_str

        with open(save_str, 'wb') as save_pickle:
            pickle.dump(fitDict, save_pickle, protocol=pickle.HIGHEST_PROTOCOL)
        
        
        
    elif coil.upper() in ['CC', 'ECE']:
        print 'Performing coil fit on', magnet.upper(), 'ECE'

        cc_param = utils.coil_datacards[magnet.upper()]['CC']
        e1_param = utils.coil_datacards[magnet.upper()]['E1']
        e2_param = utils.coil_datacards[magnet.upper()]['E2']

        cc_centre = utils.centres_dict[magnet.upper()]['mapper']['CC']
        e1_centre = utils.centres_dict[magnet.upper()]['mapper']['E1']
        e2_centre = utils.centres_dict[magnet.upper()]['mapper']['E2']

        _field = readFile(utils.coil_datacards[magnet.upper()]['CC']['30A_data'])

        if magnet.upper() == 'SSD':
            _field = fm.flip_SSD_data(_field)
        
        coilFitClass = coilfit.CoilFitClass_ECE(cc_param, e1_param, e2_param, _field, cc_centre, \
                                                e1_centre, e2_centre)
        fitDict = coilFitClass.run()
        print 'Finished with parameters: '
        for key, value in fitDict.iteritems():
            print key
            for _k, _v in value.iteritems():
                print _k, _v

        print 'Saved fit parameters at: ', save_str

        with open(save_str, 'wb') as save_pickle:
            pickle.dump(fitDict, save_pickle, protocol=pickle.HIGHEST_PROTOCOL)

    if FBfit == True:
        residField = make_resid_field(magnet.upper(), coil.upper())
        performFBfit(residField, magnet.upper(), coil.upper())
            
    return fitDict

def perform_geofit(magnet, coil, makeresid=True, save_as=None):
    if magnet.upper() not in ['SSU', 'SSD']:
        print 'Magnet unrecognised - please use SSU or SSD'
        return
    if coil.upper() not in ['M1', 'M2', 'CC', 'ECE']:
        print 'Coil unrecognised - please use M1, M2, CC or ECE'
        print '\tN.B You can not fit to the end coils individually, only to E1-CC-E2'
        return

    if coil.upper() == 'CC':
        coil = 'ECE'

    if save_as == None:
        save_str = os.path.join(utils.geofit_field_path, magnet.upper() + '_' + coil.upper() \
                                + '_geofit_default.pickle')
    else:
        save_str = os.path.join(utils.geofit_field_path, save_as)

    if coil.upper() in ['M1', 'M2']:
        print 'Performing geometrical fit on', magnet.upper(), coil.upper()
        if utils.coil_datacards[magnet.upper()][coil.upper()]['30A_data'] == None:
            print 'No data to fit to for this magnet!'
            return
        
        _centre = utils.centres_dict[magnet.upper()]['mapper'][coil.upper()]
        
        _field = readFile(utils.coil_datacards[magnet.upper()][coil.upper()]['30A_data'])
        geoFitClass = geofit.GeoFit(utils.coil_datacards[magnet.upper()][coil.upper()], \
                                              _field, _centre)
        fitDict = geoFitClass.run()
        print 'Finished with parameters: '
        for key, value in fitDict.iteritems():
            print key, value
        print 'Saved fit parameters at: ', save_str

        with open(save_str, 'wb') as save_pickle:
            pickle.dump(fitDict, save_pickle, protocol=pickle.HIGHEST_PROTOCOL)
        
        return fitDict
        
    elif coil.upper() in ['CC', 'ECE']:
        pass

    
def get_coilfit_class(magDict):
    coilFitDicts = []
    currentList = []
    _magnet = magDict['magnet']
    
    for key, item in magDict.iteritems():
        if key == 'CC':
            if item['I'] < 0.001 and item['I'] > -0.001:
                continue
            pickle_str = '%s_ECE_coilfit_default.pickle'%_magnet
            with open(os.path.join(utils.geofit_field_path, pickle_str)) as _handle:
                ece_dict = pickle.load(_handle)
                for _key, _dict in ece_dict.iteritems():
                    coilFitDicts.append(_dict)
                    currentList.append(item['I'])
        elif key in ['M1', 'M2']:
            if item['I'] < 0.001 and item['I'] > -0.001:
                continue
            pickle_str = '%s_%s_coilfit_default.pickle'%(_magnet, key)
            with open(os.path.join(utils.geofit_field_path, pickle_str)) as _handle:
                c_dict = pickle.load(_handle)
                coilFitDicts.append(c_dict)
                currentList.append(item['I'])      
          
    coilfit_class = mkfield.CalcFullField(coilFitDicts, currentList)

    return coilfit_class


def make_resid_field(magnet, coil, coilfit=True, fitDict=None, saveAs=None, _current=30.0):
    #I f*ing hate the mess that I have made this function... NEEDS CLEANING
    dataFieldStr = utils.coil_datacards[magnet.upper()][coil.upper()]['30A_data']
    
    if coil.upper() == 'CC':
        coil = 'ECE'
    if fitDict == None:
        if coilfit == True:
            fitDictStr = '%s_%s_coilfit_default.pickle'%(magnet.upper(), coil.upper())
        else:
            fitDictStr = '%s_%s_geofit_default.pickle'%(magnet.upper(), coil.upper())
    elif type(fitDict) == type('string!'):
        fitDictStr = fitDict
    elif type(fitDict) == type({}):
        fitDictStr = 'N/A'
        pass #Handle passing the actual fitDict here...

    with open(os.path.join(utils.geofit_field_path, fitDictStr), 'rb') as _file:
        fitDict = pickle.load(_file)
        
    fitDictList, currentList = [], []
    
    if coil == 'ECE':
        for key, value in fitDict.iteritems():
            fitDictList.append(value)
            currentList.append(_current)
    else:
        fitDictList.append(fitDict)
        currentList.append(_current)

    if coilfit == True:
        print 'Making residual field with coilfit using', fitDictStr, 'with data field', dataFieldStr
        calcFieldClass = mkfield.CalcFullField(fitDictList, currentList)

        dataField = readFile(dataFieldStr)

        if magnet == 'SSD':
            dataField = fm.flip_SSD_field(datafield)
        
        residualField = []
        
        for f in dataField:
            Br, Bphi, Bz = calcFieldClass.calc_full_field_at_point(f.r, f.phi, f.z)
            residualField.append(rphiz.Measurement(f.r, f.phi, f.z, f.Br - Br, f.Bphi - Bphi, \
                                                   f.Bz - Bz, f.sensorNumber))
        

    if coilfit == False:
        pass #need to implement calcgeofit class


    if saveAs == None:
        #obvs need to change this so it can handle geofit instead
        saveAs = '%s_%s_coilfit_resid.dat'%(magnet.upper(), coil.upper())

    saveAsFull = os.path.join(utils.resid_field_path, saveAs)

    fm.print_field_from_list(residualField, saveAsFull)
    
    return residualField

