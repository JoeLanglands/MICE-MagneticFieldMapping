import multiprocessing as mp
import pickle
import sys
import os

import scipy as sp
import numpy as np
import scipy.special as spec

import fbutils as _fb
from micemag.fieldmanip import polarMeasurement as rphiz
import micemag.utils as utils


#Consolidate all of this into a class to remove need for global values etc..

def getDefaultFitDict(coil, magnet):
    if coil.upper() == 'CC':
        coil = 'ECE'
        
    picklePath = os.path.join(utils.fb_pickle_path, '%s_%s_3_20_10.pickle'%(magnet, coil))
    try:
        with open(picklePath, 'rb') as _pickle:
            fitDict = pickle.load(_pickle)
        return fitDict
    except IOError:
        print 'Attempted to load pickle:', picklePath
        print 'Default FB term pickle not found! Has it been deleted?'
        print 'You need to provide the name of the pickle that you wish to use'
        sys.exit()



def applyFB_field(field, _fitDict, coil, magnet, FBonly=False, nCores=1):
    global fitDict
    fitDict = _fitDict
    global jZeros
    jZeros = _fb.genBesselZeros(fitDict['n'], fitDict['m'])
    global _FBonly
    _FBonly = FBonly
    
    _fb.centreField(field, coil, magnet)
    
    fieldPool = mp.Pool(nCores)

    field = fieldPool.map(calcFB_unpack, field)

    field.sort()
    return field


def applyFB_grid(magDict, x, y, z, Bx, By, Bz):
    global fitDict
    global jZeros
    _mag = magDict['magnet']
    _r, _phi, _z, _Br, _Bphi, _Bz = cartToPolar(x, y, z, Bx, By, Bz)
    
    if _r > 0.15:
        return Bx, By, Bz #Can't add fb terms past rmax
    
    for _coil in ['CC', 'M1', 'M2']:
        if magDict[_coil]['I'] == 0:
            continue
        if magDict[_coil]['fb'] == None:
            fitDict = getDefaultFitDict(_coil, magDict['magnet'])
        else:
            pass #need to handle non default dicts here
        jZeros = _fb.genBesselZeros(fitDict['n'], fitDict['m'])
        current_scale = magDict[_coil]['I']/30.0
        coil_centre = utils.centres_dict[_mag]['mapper'][_coil]
        Z = _z - coil_centre
        if Z < (-1.0)*fitDict['zmax'] or Z > fitDict['zmax']:
            continue
        BrFB, BphiFB, BzFB = calcBrBphiBz(_r, _phi, Z)
        _Br += BrFB*current_scale
        _Bphi += BphiFB*current_scale
        _Bz += BzFB*current_scale
        
    X, Y, Z, _Bx, _By, _Bz = polarToCart(_r, _phi, _z, _Br, _Bphi, _Bz)
    return _Bx, _By, _Bz


def cartToPolar(x, y, z, Bx, By, Bz, deg=False):
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    _matrix = np.array([[np.cos(phi), -1.0*np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    _matrix = np.transpose(_matrix)
    B = np.array([Bx, By])
    R = _matrix.dot(B)

    Br = R[0]
    Bphi = R[1]

    if deg == True:
        phi = np.degrees(phi)
        if phi < 0.0:
            phi = 360.0 + phi
        elif phi > 360:
            phi = phi - 360.0
    
    return r, phi, z, Br, Bphi, Bz


def polarToCart(r, phi, z, Br, Bphi, Bz):
    x = r*np.cos(phi)
    y = r*np.sin(phi)

    _matrix = np.array([[np.cos(phi), -1.0*np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    B = np.array([Br, Bphi])
    X = _matrix.dot(B)

    Bx = X[0]
    By = X[1]

    return x, y, z, Bx, By, Bz
    
def calcFB_unpack(point):
    global fitDict
    global _FBonly
    if point.z >= (-1.0)*fitDict['zmax'] and point.z <= fitDict['zmax'] and point.r <= 0.15:
        _Br, _Bphi, _Bz = calcBrBphiBz(point.r, np.radians(point.phi), point.z)
    else:
        _Br, _Bphi, _Bz = 0, 0, 0
    if _FBonly == False:
        return rphiz.Measurement(point.r, point.phi, point.z,\
                                 point.Br + _Br, point.Bphi+ _Bphi, point.Bz +  _Bz, \
                                 point.sensorNumber)
    elif _FBonly == True:
        return rphiz.Measurement(point.r, point.phi, point.z,\
                                 _Br, _Bphi, _Bz, point.sensorNumber)

def calcBrBphiBz(r, phi, z):
    """Calculates the Fourier Bessel field components at a point from the fitDict.
    
    phi must be in radians!
    """
    global fitDict

    Br, Bphi, Bz = 0, 0, 0
    for _n in range(fitDict['n']):
        _Brl0, _Bphil0, _Bzl0 = _calcl0terms(fitDict['A_%d_0'%_n], \
                                             fitDict['al_%d_0'%_n], \
                                             _n, r, phi, z)
        _BrE, _BphiE, _BzE = _calcEterms(fitDict['E_%d'%_n], \
                                         fitDict['ep_%d'%_n], \
                                         _n, r, phi)
        Br += _Brl0 + _BrE
        Bphi += _Bphil0 + _BphiE
        Bz += _Bzl0 #_BzE is *always* 0
        for _l in range(1, fitDict['l'] + 1):
            _BrA, _BphiA, _BzA = _calcAterms(fitDict['A_%d_%d'%(_n,_l)], \
                                             fitDict['al_%d_%d'%(_n, _l)], \
                                             _n, _l, r, phi, z)
            _BrB, _BphiB, _BzB = _calcBterms(fitDict['B_%d_%d'%(_n,_l)], \
                                             fitDict['be_%d_%d'%(_n, _l)], \
                                             _n, _l, r, phi, z)
            Br += _BrA + _BrB
            Bphi += _BphiA + _BphiB
            Bz += _BzA + _BzB
        for _m in range(1, fitDict['m'] + 1):
            _BrC, _BphiC, _BzC = _calcCterms(fitDict['C_%d_%d'%(_n, _m)], \
                                             fitDict['ga_%d_%d'%(_n, _m)], \
                                             _n, _m, r, phi, z)
            _BrD, _BphiD, _BzD = _calcDterms(fitDict['D_%d_%d'%(_n, _m)], \
                                             fitDict['de_%d_%d'%(_n, _m)], \
                                             _n, _m, r, phi, z)
            Br += _BrC + _BrD
            Bphi += _BphiC + _BphiD
            Bz += _BzC + _BzD
        
    return Br, Bphi, Bz

    
def _calcAterms(A, al, n, l, r, phi, z):
    global fitDict
    sc = np.pi/(fitDict['zmax']*1.2)
    Br = A*spec.ivp(n, l*sc*r)*np.cos(n*phi + al)*np.sin(l*sc*z)
    if r == 0:
        Bphi = 0
    else:
        Bphi = (-1)*A*(1/(l*sc*r))*spec.iv(n, l*sc*r)*np.sin(n*phi + al)*np.sin(l*sc*z)
    Bz = A*spec.iv(n, l*sc*r)*np.cos(n*phi + al)*np.cos(l*sc*z)
    return Br, Bphi, Bz

def _calcBterms(B, be, n, l, r, phi, z):
    global fitDict
    sc = np.pi/(fitDict['zmax']*1.2)
    Br = B*spec.ivp(n, l*sc*r)*np.cos(n*phi + be)*np.cos(l*sc*z)
    if r == 0:
        Bphi = 0
    else:
        Bphi = (-1)*B*(1/(l*sc*r))*spec.iv(n, l*sc*r)*np.sin(n*phi + be)*np.cos(l*sc*z)
    Bz =  (-1)*B*spec.iv(n, l*sc*r)*np.cos(n*phi + be)*np.sin(l*sc*z)
    return Br, Bphi, Bz

def _calcl0terms(A, al, n, r, phi, z):
    if r == 0 and n == 0:
        Br, Bphi = 0, 0
    else:
        Br = A*n*np.power(r, n-1)*np.cos(n*phi + al)*z
        Bphi = (-1)*A*np.power(r, n-1)*np.sin(n*phi + al)*z
    Bz = A*np.power(r, n)*np.cos(n*phi + al)
    return Br, Bphi, Bz

def _calcCterms(C, ga, n, m, r, phi, z):
    global jZeros
    global fitDict
    sc = jZeros[n][m]/fitDict['rmax']
    Br = C*spec.jvp(n, sc*r)*np.cos(n*phi + ga)*np.sinh(sc*z)
    if r == 0:
        Bphi = 0
    else:
        Bphi = (-1)*C*(1/(sc*r))*spec.jv(n, sc*r)*np.sin(n*phi + ga)*np.sinh(sc*z)
    Bz = C*spec.jv(n, sc*r)*np.cos(n*phi + ga)*np.cosh(sc*z)
    return Br, Bphi, Bz

def _calcDterms(D, de, n, m, r, phi, z):
    global jZeros
    global fitDict
    sc = jZeros[n][m]/fitDict['rmax']
    Br = D*spec.jvp(n, sc*r)*np.cos(n*phi + de)*np.cosh(sc*z)
    if r == 0:
        Bphi = 0
    else:
        Bphi = (-1)*D*(1/(sc*r))*spec.jv(n, sc*r)*np.sin(n*phi + de)*np.cosh(sc*z)
    Bz = D*spec.jv(n, sc*r)*np.cos(n*phi + de)*np.sinh(sc*z)
    return Br, Bphi, Bz

def _calcEterms(E, ep, n, r, phi):
    if r == 0 and n == 0:
        Br, Bphi = 0, 0
    else:
        Br = E*n*np.power(r, n-1)*np.cos(n*phi + ep)
        Bphi = (-1)*E*n*np.power(r, n-1)*np.sin(n*phi + ep)
    Bz = 0
    return Br, Bphi, Bz
