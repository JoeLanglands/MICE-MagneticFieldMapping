import sys

import scipy as sp
import numpy as np
import scipy.special as spec

import micemag.utils as utils
from micemag.fieldmanip.fieldManipulation import shiftField


def centreField(field, coil, magnet, undo=False):
    """Function that smartly shifts the data so that z=0 is where the coil centre is."""
    field.sort()

    if field[0].z > 8000:
        data_set = 'survey'
    else:
        data_set = 'mapper'

    if magnet in ['ssu', 'SSU', 'upstream']:
        shift = utils.ssu_centre_dict[data_set][coil]
    elif magnet in ['ssd', 'SSD', 'downstream']:
        shift = utils.ssd_centre_dict[data_set][coil]
        
    if undo == False:
        mult = -1.0
    elif undo == True:
        mult = 1.0
        
    shiftField(field, mult*shift)

    
def mergeDicts(*_dicts):
    """Function that takes any multiple of dictionaries as arguments and merges them together."""
    result = {}
    for d in _dicts:
        result.update(d)
    return result

def genBesselZeros(n, m):
    """Function that generates m roots of each order of Bessel function.

    Args:
        n: Number of orders of Bessel functions.
        m: Number of desired roots of each order.

    Returns:
        A 2D array filled with the roots of each order that can be accessed easily.
        The mth zero of the nth Bessel function can be accessed via besselZeros[n][m].
        It also fills the arrays with a leading '0' element so that [m-1] IS NOT used
        as an attempt to make it more user friendly.
    """
    besselZeros = []
    for _n in range(0, n+1):
        besselZeros.append(np.append([0], spec.jn_zeros(_n, m)))
    return besselZeros


def _chunkIter(seq, chSize=2):
    """Helper function for iterating through the generated arguments in pairs."""
    return (seq[pos: pos + chSize] for pos in xrange(0, len(seq), chSize))


class _FourierFcnArgGen:
    def __init__(self, n, l):
        """Class to trick minuit into taking a function with a variable number of parameters.

        This class is for the Fourier terms of the FB expansion.  It is designed to work with
        the 'FourierFitFcn' class which actually does the work.
        """
        self.co_varnames = ()
        for _n in range(n):
            for _l in range(1, l+1):
                self.co_varnames += ('A_%d_%d'%(_n, _l), \
                                     'al_%d_%d'%(_n, _l), \
                                     'B_%d_%d'%(_n,_l), \
                                     'be_%d_%d'%(_n, _l))
        for _n in range(n):
            self.co_varnames += ('A_%d_0'%_n, 'al_%d_0'%_n)

        self.co_argcount = len(self.co_varnames)

class FourierFitFcn:
    def __init__(self, data, n, l, zmax=1.0, rmax=0.15, sigma=0.001, verbose=True, **kwargs):
        #Format of data is list of tuples
        #e.g [(r, phi, z, Bz), (r, phi, z, Bz), ...]
        #This will be sorted in the full fit class
        self.func_code = _FourierFcnArgGen(n, l)
        self.data = data
        self.n = n
        self.l = l
        self.zmax = zmax
        self.rmax = rmax
        self.chiSq = 0
        self.sigma = sigma
        self.Coeffs = None
        self.verbose = verbose
        self.DoF = len(data) - self.func_code.co_argcount #Degrees of freedom

    def __call__(self, *args):
        self.chiSq = 0
        self.Coeffs = args # save the args so we can calculate the terms with the fitted params

        for _data in self.data:
            _r = _data[0]
            _phi = _data[1]
            _z = _data[2]
            _dataBz = _data[3]
            
            modelBz = self.calcFourierTerms(_r, _phi, _z)
            self._calcSqRes(_dataBz, modelBz)
    
        if self.verbose == True:
            print self.chiSq, self.DoF, self.chiSq/self.DoF
        return self.chiSq
            
   
    def _calcAterms(self, A, al, n, l, r, phi, z):
        sc = np.pi/(self.zmax*1.2)
        return A*spec.iv(n, l*sc*r)*np.cos(n*phi + al)*np.cos(l*sc*z)

    def _calcBterms(self, B, be, n, l, r, phi, z):
        sc = np.pi/(self.zmax*1.2)
        return (-1)*B*spec.iv(n, l*sc*r)*np.cos(n*phi + be)*np.sin(l*sc*z)

    def _calcl0terms(self, A, al, n, r, phi, z):
        return A*np.power(r, n)*np.cos(n*phi + al)

    def _calcAtermsBr(self, A, al, n, l, r, phi, z):
        sc = np.pi/(self.zmax*1.2)
        return A*spec.ivp(n, l*sc*r)*np.cos(n*phi + al)*np.sin(l*sc*z)

    def _calcBtermsBr(self, B, be, n, l, r, phi, z):
        sc = np.pi/(self.zmax*1.2)
        return B*spec.ivp(n, l*sc*r)*np.cos(n*phi + be)*np.cos(l*sc*z)

    def _calcl0termsBr(self, A, al, n, r, phi, z):
        return A*n*np.power(r, n-1)*np.cos(n*phi + al)*z

    def _calcSqRes(self, dataBz, modelBz):
        sqRes = (modelBz - dataBz)**2/(self.sigma**2)
        self.chiSq += sqRes

    def calcFourierTerms(self, r, phi, z, comp='Bz'):
        _sum = 0
        for i, (_arg, _argName) in enumerate(zip(_chunkIter(self.Coeffs), \
                                                 _chunkIter(self.func_code.co_varnames))):
            tmpStr = _argName[0].split('_')
            _n = int(tmpStr[1])
            _l = int(tmpStr[2])
        
            if i%2 == 0 and i < (self.l*2*self.n):
                if comp == 'Bz':
                    _sum += self._calcAterms(_arg[0], _arg[1], _n, _l, r, phi, z)
                elif comp == 'Br':
                    _sum += self._calcAtermsBr(_arg[0], _arg[1], _n, _l, r, phi, z)
            elif i%2 == 1 and i < (self.l*2*self.n):
                if comp == 'Bz':
                    _sum += self._calcBterms(_arg[0], _arg[1], _n, _l, r, phi, z)
                elif comp == 'Br':
                    _sum += self._calcBtermsBr(_arg[0], _arg[1], _n, _l, r, phi, z)
            elif i >= (self.l*2*self.n):
                if comp == 'Bz':
                    _sum += self._calcl0terms(_arg[0], _arg[1], _n, r, phi, z)
                elif comp == 'Br':
                    _sum += self._calcl0termsBr(_arg[0], _arg[1], _n, r, phi, z)

        return _sum

    def setData(self, data):
        self.data == data
   


class _HyperbolicFcnArgGen:
    def __init__(self, n, m):
        """Class to trick minuit into taking a function with a variable number of parameters.

        This class is for finding the Hyperbolic terms of the FB expansion.  It is designed to
        work with the 'HyperbolicFitFcn' class which actually does the work.
        """
        self.co_varnames = ()
        for _n in range(n):
            for _m in range(1, m+1):
                self.co_varnames += ('C_%d_%d'%(_n, _m), \
                                     'ga_%d_%d'%(_n, _m), \
                                     'D_%d_%d'%(_n,_m), \
                                     'de_%d_%d'%(_n, _m))

        self.co_argcount = len(self.co_varnames)


class HyperbolicFitFcn:
    def __init__(self, data, n, m, zmax=1.0, rmax=0.15, sigma=0.001, verbose=True, **kwargs):
        self.func_code = _HyperbolicFcnArgGen(n, m)
        self.data = data
        self.n = n
        self.m = m
        self.zmax = zmax
        self.rmax = rmax
        self.chiSq = 0
        self.sigma = sigma
        self.Coeffs = None
        self.jZeros = genBesselZeros(n, m)
        self.verbose = verbose
        self.DoF = None
        

    def __call__(self, *args):
        if self.data == None:
            print 'There is no data to fit to yet!!'
            return None
        self.chiSq = 0
        self.Coeffs = args # save the args so we can calculate the terms with the fitted params

        for _data in self.data:
            _r = _data[0]
            _phi = _data[1]
            _z = _data[2]
            _dataBz = _data[3]
            
            modelBz = self.calcHypTerms(_r, _phi, _z)
            self._calcSqRes(_dataBz, modelBz)

        if self.verbose == True:
            print self.chiSq, self.DoF, self.chiSq/self.DoF
        return self.chiSq

    def _calcCterms(self, C, ga, n, m, r, phi, z):
        sc = self.jZeros[n][m]/self.rmax
        return C*spec.jv(n, sc*r)*np.cos(n*phi + ga)*np.cosh(sc*z)

    def _calcDterms(self, D, de, n, m, r, phi, z):
        sc = self.jZeros[n][m]/self.rmax
        return D*spec.jv(n, sc*r)*np.cos(n*phi + de)*np.sinh(sc*z)

    def _calcCtermsBr(self, C, ga, n, m, r, phi, z):
        sc = self.jZeros[n][m]/self.rmax
        return C*spec.jvp(n, sc*r)*np.cos(n*phi + ga)*np.sinh(sc*z)

    def _calcDtermsBr(self, D, de, n, m, r, phi, z):
        sc = self.jZeros[n][m]/self.rmax
        return D*spec.jvp(n, sc*r)*np.cos(n*phi + de)*np.cosh(sc*z)

    def _calcSqRes(self, dataBz, modelBz):
        sqRes = (modelBz - dataBz)**2/(self.sigma**2)
        self.chiSq += sqRes

    def calcHypTerms(self, r, phi, z, comp='Bz'):
        _sum = 0
        for i, (_arg, _argName) in enumerate(zip(_chunkIter(self.Coeffs), \
                                                 _chunkIter(self.func_code.co_varnames))):
            tmpStr = _argName[0].split('_')
            _n = int(tmpStr[1])
            _m = int(tmpStr[2])
        
            if i%2 == 0:
                if comp == 'Bz':
                    _sum += self._calcCterms(_arg[0], _arg[1], _n, _m, r, phi, z)
                elif comp == 'Br':
                    _sum += self._calcCtermsBr(_arg[0], _arg[1], _n, _m, r, phi, z)
            elif i%2 == 1:
                if comp == 'Bz':
                    _sum += self._calcDterms(_arg[0], _arg[1], _n, _m, r, phi, z)
                elif comp == 'Br':
                    _sum += self._calcDtermsBr(_arg[0], _arg[1], _n, _m, r, phi, z)

        return _sum

    def setData(self, data):
        #Need function to set data so we can initialize this class with data = None.
        #Then it can be set with this function after it has been found from the fourier
        #terms and the field data. It gives a little more flexibility.
        self.DoF = len(data) - self.func_code.co_argcount #Degrees of freedom
        self.data = data


class _MultipoleFcnArgGen:
    def __init__(self, n):
        """Class to trick minuit into taking a function with a variable number of parameters.

        This class is for finding the Multipole terms of the FB expansion.  It is designed to
        work with the 'MultipoleFitFcn' class which actually does the work.
        """
        self.co_varnames = ()
        for _n in range(n):
            self.co_varnames += ('E_%d'%_n, \
                                 'ep_%d'%_n)

        self.co_argcount = len(self.co_varnames)

class MultipoleFitFcn:
    def __init__(self, data, n, rmax=0.15, sigma=0.001, verbose=True, **kwargs):
        self.func_code = _MultipoleFcnArgGen(n)
        self.n = n
        self.sigma = sigma
        self.Coeffs = None
        self.data = data #in a similar format but with [(phi, avBr),...]
                         #no need for z or r!
        self.rmax = rmax
        self.verbose = verbose
        self.DoF = None

    def __call__(self, *args):
        if self.data == None:
            print 'There is no data to fit to yet!!'
            return None
        self.chiSq = 0
        self.Coeffs = args # save the args so we can calculate the terms with the fitted params

        for _data in self.data:
            _phi = _data[0]
            _dataBr = _data[1]
            
            modelBr = self.calcMultipoleTerms(self.rmax, _phi)
            self._calcSqRes(_dataBr, modelBr)
        if self.verbose == True:
            print self.chiSq, self.DoF, self.chiSq/self.DoF
        return self.chiSq


    def calcMultipoleTerms(self, r, phi): #Doesn't depend on z at all
        _sum = 0
        for i, (_arg, _argName) in enumerate(zip(_chunkIter(self.Coeffs), \
                                                 _chunkIter(self.func_code.co_varnames))):
            tmpStr = _argName[0].split('_')
            _n = int(tmpStr[1])
            
            _sum += self._calcEterms(_arg[0], _arg[1], _n, r, phi)

        return _sum

    def _calcEterms(self, E, ep, n, r, phi):
        try:
            return E*n*np.power(r, n-1)*np.cos(n*phi + ep)
        except RuntimeWarning:
            return 0
    def _calcSqRes(self, dataBr, modelBr):
        sqRes = (modelBr - dataBr)**2/(self.sigma**2)
        self.chiSq += sqRes

    def setData(self, data):
        self.DoF = len(data) - self.func_code.co_argcount #Degrees of freedom
        self.data = data


if __name__=='__main__':
    centreField(0, 0 ,0)
