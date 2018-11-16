import pickle
import sys
import os

import scipy as sp
import numpy as np
import iminuit as minuit

import fbutils as _fb
import micemag.utils as _paths


class FBfitClass:
    def __init__(self, field, coil, magnet, zmax=1.8, rmax=0.15, n=2, l=20, m=10, \
                 verbose=True, saveAs=None):
        self.field = field
        _fb.centreField(self.field, coil, magnet)
        self.n = n
        self.l = l
        self.m = m
        self.coil = coil
        self.magnet = magnet
        self.fitDict = {'n': n, 'l': l, 'm': m, 'zmax': zmax, 'rmax': rmax}
        self.zmax = zmax
        self.rmax = rmax
        self.fouData, self.zstart, self.zend = self._getFourierData()
        self.FourFitClass = _fb.FourierFitFcn(self.fouData, n, l, zmax=zmax, \
                                              rmax=rmax,\
                                              verbose=verbose)
        self.HypFitClass = _fb.HyperbolicFitFcn(None, n, m, zmax=zmax, \
                                                rmax=rmax, \
                                                verbose=verbose)
        self.MultFitClass = _fb.MultipoleFitFcn(None, n, rmax=rmax, verbose=verbose)
        self.verbose = verbose
        if saveAs == None:
            self.saveAs = magnet + '_' + coil + '_%d_%d_%d.pickle'%(n, l, m) 
        else:
            self.saveAs = saveAs

    def run(self, save=True): #The main function that runs everything
        self._runFourierFit()

        self._getHyperbolicData()
        self._runHyperbolicFit()

        self._getMultipoleData()
        self._runMultipoleFit()

        self.saveDict(self.saveAs)
        
        return self.fitDict

    def saveDict(self, saveAs):
        if saveAs[-7:] != '.pickle':
            saveAs += '.pickle'
        with open(os.path.join(_paths.fb_pickle_path, saveAs), 'wb') as _pickle:
            pickle.dump(self.fitDict, _pickle, protocol=pickle.HIGHEST_PROTOCOL)
                                             
    def _getFourierData(self):
        _fourierData = []
        for point in self.field:
            if point.z >= -1*self.zmax and point.z <= self.zmax and point.r == self.rmax:
               _fourierData.append((point.r, np.radians(point.phi),\
                                     point.z, point.Bz))
        _z = []
        for d in _fourierData:
            _z.append(d[2])
        zstart = min(_z)
        zend = max(_z)
        return _fourierData, zstart, zend

    def _runFourierFit(self):
        if self.verbose == True:
            _min = minuit.Minuit(self.FourFitClass)
        else:
            _min = minuit.Minuit(self.FourFitClass, pedantic=False, print_level=0)
        _min.migrad()

        self.fitDict = _fb.mergeDicts(self.fitDict, _min.values)

    def _getHyperbolicData(self):
        _hyperbolicData = []
        for i in self.field:
            if i.z > self.zstart - 0.001 and i.z < self.zstart + 0.001:
                fourBz = self.FourFitClass.calcFourierTerms(i.r, np.radians(i.phi), i.z)
                _hyperbolicData.append((i.r, np.radians(i.phi), i.z, \
                                        i.Bz - fourBz))
            if i.z > self.zend - 0.001 and i.z < self.zend + 0.001:
                fourBz = self.FourFitClass.calcFourierTerms(i.r, np.radians(i.phi), i.z)
                _hyperbolicData.append((i.r, np.radians(i.phi), i.z, \
                                        i.Bz - fourBz))
        self.HypFitClass.setData(_hyperbolicData)
        return _hyperbolicData

    def _runHyperbolicFit(self):
        #This MUST be called AFTER _runFourierFit() and AFTER _getHyperbolicData
        if self.verbose == True:
            _min = minuit.Minuit(self.HypFitClass)
        else:
            _min = minuit.Minuit(self.HypFitClass, pedantic=False, print_level=0)
        _min.migrad()

        self.fitDict = _fb.mergeDicts(self.fitDict, _min.values)
    

    def _getMultipoleData(self):
        multipoleData = []
        uniquePhis = {} # {phi: sumBr, phi2: sumBr2, ...}
        for point in self.field:
            if point.r == self.rmax and point.z >= self.zstart and point.z <= self.zend:
                _r = point.r
                _phi = np.radians(point.phi)
                _z = point.z
                if point.phi in uniquePhis:
                    uniquePhis[point.phi] += point.Br \
                                             - self.FourFitClass.calcFourierTerms(_r, _phi, _z, \
                                                                                  comp='Br') \
                                             - self.HypFitClass.calcHypTerms(_r, _phi, _z, \
                                                                             comp='Br')
                else:
                    uniquePhis[point.phi] = point.Br \
                                            - self.FourFitClass.calcFourierTerms(_r, _phi, _z, \
                                                                                 comp='Br') \
                                            - self.HypFitClass.calcHypTerms(_r, _phi, _z, \
                                                                            comp='Br')
                    
        for key, value in uniquePhis.iteritems():
            multipoleData.append((key, np.mean(value)))

        self.MultFitClass.setData(multipoleData)
                                 
        
    def _runMultipoleFit(self):
        if self.verbose == True:
            _min = minuit.Minuit(self.MultFitClass)
        else:
            _min = minuit.Minuit(self.MultFitClass, pedantic=False, print_level=0)
        _min.migrad()

        self.fitDict = _fb.mergeDicts(self.fitDict, _min.values)

    def getFitDict(self):
        return self.fitDict
