import sys

import matplotlib.pyplot as pyplot
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import scipy.interpolate as sp

from micemag.fieldmanip import polarMeasurement as rphiz

np.set_printoptions(threshold=np.nan)

def wireFrame(field, field2=None, comp='Bz', unit='T', _R=(0.15, 0.15), zlim=None, labels=('field1', 'field2'), **kwargs):
    """
    This will be a generic version of ContourRPhiplot.py so that these plots can be made easily

    """
    
    if unit == 'T':
        _mult = 1
    elif unit == 'mT':
        _mult = 1000
    elif unit == 'G':
        _mult = 10000
    else: #if it is ballsed up then just leave it as Tesla
        unit = 'T'
        _mult = 1
        
    field =  roundField(field, Z=True)

    field.sort()

    if zlim == None:
        zlim = (field[0].z, field[-1].z)
    
  
    if 'rphi' in kwargs:
        X, Y, Z = _shapeDataRPHI(field, comp=comp)
    else:
        X, Y, Z = _shapeData(field, _R[0], comp, _mult, zlim)
   
    
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
        
    if comp == 'Bphi':
        z_label = r'$B_{\phi}$ [%s]'%unit
    elif comp == 'B':
        z_label = r'|$B$| [%s]'%unit
    else:
        z_label = r'$\Delta B_{%s}$'%comp[1] + ' [%s]'%unit
        
    ax.set_xlabel('z [m]', fontsize=14)
    ax.set_ylabel(r'$\phi$ [deg]',fontsize=14)
    ax.set_zlabel('        ' +z_label, fontsize=14)
    if 'title' in kwargs:
        ax.set_title(kwargs['title'], fontsize=14)

    ax.plot_wireframe(X, Y, Z, color='m', alpha=0.5, label=labels[0])
    

                
    if field2 != None:
        field2 = roundField(field2, Z=True) 

        X2, Y2, Z2 = _shapeData(field2, _R[1], comp, _mult, zlim)
            
        ax.plot_wireframe(X2, Y2, Z2, color='r', alpha=0.5, label=labels[1])
        
        
    ax.zaxis.set_rotate_label(False)
    ax.zaxis._axinfo['label']['space_factor'] = 10
    pyplot.tight_layout()
    ax.view_init(30,-75)
    if 'saveAs' in kwargs:
        pyplot.savefig(kwargs['saveAs'])

    ax.legend()    
    
    pyplot.show()


def _shapeData(field, _R=0.15, comp='Bz', _mult=1, zlim=(-1.0, 1.0)):
    field.sort()
    #checkGrid(field)
    zStart = np.around(field[0].z, 2)
    PhiList = []
    _phi_skip = None
    for f in field:
        if np.around(f.r, 2) == _R and np.around(f.phi,0) in PhiList:
            _phi_skip = np.around(f.phi,0)
            break
        if np.around(f.r, 2) == _R:
            PhiList.append(np.around(f.phi,0))
        if np.around(f.z, 2) != zStart:
            break

    _z, _phi, _B = [], [], []
    skipPhi = False
    for f in field:
        if np.around(f.phi,0) == _phi_skip and skipPhi == True:
            skipPhi = False
            continue
        if np.around(f.r,2) == _R and np.around(f.z,2) > zlim[0] and np.around(f.z,2) < zlim[1]:
            #print f.z, f.phi
            _z.append(np.around(f.z, 2))
            _phi.append(np.around(f.phi,0))
            if comp == 'Bz':
                _B.append(f.Bz*_mult)
            elif comp == 'Br':
                _B.append(f.Br*_mult)
            elif comp == 'Bphi':
                _B.append(f.Bphi*_mult)
            if f.phi == _phi_skip:
                skipPhi = True
    _z.sort()
    
    _z = np.array(_z)
    _phi = np.array(_phi)
    _B = np.array(_B)
    
    cols = np.unique(_phi).shape[0]

    _X = _z.reshape(-1, cols)
    _Y = _phi.reshape(-1, cols)
    _Z = _B.reshape(-1, cols)

    
    return _X, _Y, _Z


def _shapeDataRPHI(field, comp='Bz', _mult=1):
    field.sort()
  
    _z, _phi, _B = [], [], []
    skipPhi = False
    for f in field:
        _z.append(np.around(f.r,2))
        _phi.append(np.around(f.phi, 0))
        if comp == 'Bz':
            _B.append(f.Bz*_mult)
        elif comp == 'Br':
            _B.append(f.Br*_mult)
        elif comp == 'Bphi':
            _B.append(f.Bphi*_mult)
        elif comp == 'B':
            _B.append(f.B*_mult)
        

    _z.sort()

    print len(_z), len(_phi), len(_B)
    
    _z = np.array(_z)
    _phi = np.array(_phi)
    _B = np.array(_B)
    
    cols = np.unique(_phi).shape[0]
   
    X = _z.reshape(-1, cols)
    Y = _phi.reshape(-1, cols)
    Z = _B.reshape(-1, cols)

    return X, Y, Z



def roundField(field, zAround=3, **kwargs):
    """
    This function basically truncates and rounds number in the data to get rid of the 0.000012345 from stuff.
    When rotated a point is sometimes r = 0.11999 instead of 0.12 which is annoying.
    Also it rounds the z value intelligently so instead of z being 0.1000009876 it is just 0.10.
    """
    resultField = []
    for f in field:
        if f.r < 0.0001:
            PHI = 0.0
        else:
            PHI = np.around(f.phi, 0)
        if 'Z' in kwargs:
            if kwargs['Z'] == True:
                Z = np.around(f.z, zAround)
            elif kwargs['Z'] == False:
                Z = f.z
        else:
            Z = f.z
        resultField.append(rphiz.Measurement(np.around(f.r, 2), PHI, Z, f.Br, f.Bphi, f.Bz, f.sensorNumber))

    return resultField
