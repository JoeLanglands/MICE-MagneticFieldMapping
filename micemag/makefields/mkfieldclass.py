import copy_reg
import types
import sys
import multiprocessing as mp

import numpy as np

import micemag.fieldmanip.polarMeasurement as rphiz
try:
    import makefield_cpp as mkcpp #C++ lib
except ImportError as err:
    print err.__class__.__name__ + ': ' + err.message
    print 'Has the cpp file been compiled? Try running make in the top level directory'
    sys.exit()
    
#Below function from: https://stackoverflow.com/questions/27318290/why-can-i-pass-an-instance-method-to-multiprocessing-process-but-not-a-multipro
#Allows class methods to be pickled so that they can be used by multiprocessing pool
def _reduce_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)
copy_reg.pickle(types.MethodType, _reduce_method)


class CalcField(object):
    def __init__(self, length=0, rInner=0, rOuter=0, nLayers=0, nTurns=0, current=0, centrePos=0,
                 thetaX=0, thetaY=0, px=0, py=0, fitParamDict=None):
        #static typing as they're gonna be fed into c++
        #These things are physical properties of the coil
        if fitParamDict == None:
            self.length = float(length)
            self.rInner = float(rInner)
            self.rOuter = float(rOuter)
            self.nLayers = int(nLayers)
            self.nTurns = int(nTurns)

            #These things can change
            self.current = float(current)
            self.centre = float(centrePos)
            self.thetaX = float(thetaX)
            self.thetaY = float(thetaY)
            self.px = float(px)
            self.py = float(py)

        else:
            self.set_with_dict(fitParamDict)
            self.current = current

    def set_with_dict(self, fitParamDict): #So the coilfit parameter dicts can be used
        try:
            self.length = fitParamDict['length']
            self.rInner = fitParamDict['rInner']
            self.rOuter = fitParamDict['rOuter']
            self.nLayers = int(fitParamDict['nLayers'])
            self.nTurns = int(fitParamDict['nTurns'])
            self.centre = fitParamDict['centre']
            self.thetaX = fitParamDict['thetaX']
            self.thetaY = fitParamDict['thetaY']
            self.px = fitParamDict['px']
            self.py = fitParamDict['py']
        except KeyError: #needs actual exception handling written...
            return
        
            
    def set_current(self, current):
        self.current = float(current)

    def set_centre_pos(self, centrePos):
        self.centre = float(centrePos)

    def set_angles_offsets(self, thetaX, thetaY, px, py):
        self.thetaX = float(thetaX)
        self.thetaY = float(thetaY)
        self.px = float(px)
        self.py = float(py)

    def calc_field_at_point_xyz(self, x, y, z):
        Bx, By, Bz = mkcpp.get_field_at_point_xyz(self.current, self.centre, self.length, \
                                                  self.rInner, self.rOuter, self.nLayers, \
                                                  self.nTurns, self.thetaX, self.thetaY, \
                                                  self.px, self.py, x, y, z)
        return Bx, By, Bz
    
    def calc_field_at_point(self, r, phi, z): #PHI IN DEGREES (converted here)
        Br, Bphi,  Bz = mkcpp.get_field_at_point(self.current, self.centre, self.length, \
                                                 self.rInner, self.rOuter, self.nLayers, \
                                                 self.nTurns, self.thetaX, self.thetaY, \
                                                 self.px, self.py, r, np.radians(phi), z)
        return Br, Bphi, Bz
            
    def calc_field_measurement(self, point): #Phi is in degrees cos it's from a measurement class
        Br, Bhi, Bz = mkcpp.get_field_at_point(self.current, self.centre, self.length, \
                                               self.rInner, self.rOuter, self.nLayers, \
                                               self.nTurns, self.thetaX, self.thetaY, \
                                               self.px, self.py, point.r, \
                                               np.radians(point.phi),  point.z)

        result_point = rphiz.Measurement(point.r, point.phi, point.z, Br, point.Bphi, Bz)
        return result_point

   
class CalcFullField(object):
    """
    Convenience class that uses multiple CalcField classes to build the field from mulitple coils.
    Easy to initialize from a list of the coilfit dictionaries.
    """
    def __init__(self, coilFitDicts, currentList):
        if len(coilFitDicts) != len(currentList):
            print 'You have set', len(currentList), 'currents for', len(coilFitDicts), 'coils!'
            sys.exit()
        self.calc_classes = []
        
        for n, _dict in enumerate(coilFitDicts):
            self.calc_classes.append(CalcField())
            self.calc_classes[n].set_with_dict(_dict)
            self.calc_classes[n].set_current(currentList[n])


    def calc_full_field_at_point(self, r, phi, z):
        sumBr, sumBphi, sumBz = 0, 0, 0
        for cls in self.calc_classes:
            tmpbr, tmpbp, tmpbz = cls.calc_field_at_point(r, phi, z)
            sumBr += tmpbr
            sumBphi += tmpbp
            sumBz += tmpbz

        return sumBr, sumBphi, sumBz
            
    def calc_full_field_at_point_xyz(self, x, y, z):
        sumBx, sumBy, sumBz = 0, 0, 0
        for cls in self.calc_classes:
            tmpbx, tmpby, tmpbz = cls.calc_field_at_point_xyz(x, y, z)
            sumBx += tmpbx
            sumBy += tmpby
            sumBz += tmpbz

        return sumBx, sumBy, sumBz
