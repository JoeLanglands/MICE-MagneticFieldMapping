import sys
import os
import pickle

import scipy as sp
import numpy as np
import iminuit as minuit

import micemag.makefields.mkfieldclass as mkfield
import geofit_utils as gu

    
class CoilFitClass(object):
    def __init__(self, magParamDict, field=None, centre_pos=None, sigma=0.002):
        #These things will not change
        self.data_field = gu.cut_field_for_fit(field)
        if field != None:
            self.current = 30.0*gu.get_polarity(field)
            if field[0].ID == 'Cartesian Data':
                print 'Need to pass field as polar coordinates to ', self.__class__.__name__
        self.nTurns = magParamDict['nTurns']
        self.nLayers = magParamDict['nLayers']
        self.sigma = sigma

        #These things could change in the fit and are used as seeds
        self.length = magParamDict['length']
        self.rInner = magParamDict['rInner']
        self.rOuter = magParamDict['rOuter']
        if centre_pos != None:
            self.centre = centre_pos #If you want to run a fit a seed centre_pos should be given
        
        try:
            self.thetaX = magParamDict['thetaX']
            self.thetaY = magParamDict['thetaY']
            self.px = magParamDict['px']
            self.py = magParamDict['py']
            self.centre = magParamDict['centre'] 
        except KeyError:
            #If they do not exist then set them as 0 for seeds for the fit
            self.thetaX = 0
            self.thetaY = 0
            self.px = 0
            self.py = 0
            
            
    def set_fit_data(self, field):
        self.data_field = gu.cut_field_for_fit(field)
        
    def run(self):
        if self.data_field == None:
            print 'No data to fit too!'
            return 
        
        print 'Fitting to', len(self.data_field), 'data points'
        #set sensible limits for the coil parameters otherwise it can cock up completely
        rIlow = 0.256 #Simply from the bore
        rIhi = self.rInner*1.15 #15% increase
        rOlow = self.rOuter*0.85 #15% decrease
        rOhi = self.rOuter*1.15
        lenLow = self.length*0.85
        lenHi = self.length*1.15

        minimizer = minuit.Minuit(self._min_func, rInner=self.rInner, limit_rInner=(rIlow, rIhi), \
                                  rOuter=self.rOuter, limit_rOuter=(rOlow, rOhi), \
                                  length=self.length, limit_length=(lenLow, lenHi), \
                                  centre=self.centre, thetaX=self.thetaX, thetaY=self.thetaY, \
                                  px=self.px, py=self.py, pedantic=False)

        minimizer.migrad()

        fitDict = minimizer.values
        fitDict['nLayers'] = self.nLayers
        fitDict['nTurns'] = self.nTurns

        return fitDict
            

    def _min_func(self, rInner, rOuter, length, centre, thetaX, thetaY, px, py):
        calc_class = mkfield.CalcField(length, rInner, rOuter, self.nLayers, self.nTurns, \
                                       self.current, centre, thetaX, thetaY, px, py)

        sum_br_res = 0
        sum_bp_res = 0
        sum_bz_res = 0
        
        for point in self.data_field:
            tmp_point = calc_class.calc_field_measurement(point)
            sum_br_res += (point.Br - tmp_point.Br)**2/(self.sigma*self.sigma)
            sum_bp_res += (point.Bphi - tmp_point.Bphi)**2/(self.sigma*self.sigma)
            sum_bz_res += (point.Bz - tmp_point.Bz)**2/(self.sigma*self.sigma)

        chi_sq = sum_br_res + sum_bp_res + sum_bz_res
        print chi_sq/(len(self.data_field)-8)
        return chi_sq

   
        

class CoilFitClass_ECE(object):
    def __init__(self, CCParamDict, E1ParamDict, E2ParamDict, field=None, CCcentre=0, E1centre=0,  \
                 E2centre=0, sigma=0.002):
        self.data_field = gu.cut_field_for_fit(field)
        self.minuit_kwargs = {}

        self.nLayersDict = {}
        self.nTurnsDict = {}

        self.setup_CC(CCParamDict)
        self.setup_E1(E1ParamDict)
        self.setup_E2(E2ParamDict)

        self.CC_centre = CCcentre
        self.E1_centre = E1centre
        self.E2_centre = E2centre

        self.sigma = sigma
   
    def setup_CC(self, paramDict):
        self.CC_nLayers = paramDict['nLayers']
        self.CC_nTurns = paramDict['nTurns']
        self.CC_length = paramDict['length']
        self.CC_rInner = paramDict['rInner']
        self.CC_rOuter = paramDict['rOuter']
        self.nLayersDict['CC'] = paramDict['nLayers']
        self.nTurnsDict['CC'] = paramDict['nTurns']
        
    def setup_E1(self, paramDict):
        self.E1_nLayers = paramDict['nLayers']
        self.E1_nTurns = paramDict['nTurns']
        self.E1_length = paramDict['length']
        self.E1_rInner = paramDict['rInner']
        self.E1_rOuter = paramDict['rOuter']
        self.nLayersDict['E1'] = paramDict['nLayers']
        self.nTurnsDict['E1'] = paramDict['nTurns']
        
    def setup_E2(self, paramDict):
        self.E2_nLayers = paramDict['nLayers']
        self.E2_nTurns = paramDict['nTurns']
        self.E2_length = paramDict['length']
        self.E2_rInner = paramDict['rInner']
        self.E2_rOuter = paramDict['rOuter']
        self.nLayersDict['E2'] = paramDict['nLayers']
        self.nTurnsDict['E2'] = paramDict['nTurns']
        
    def set_E1_centre(self, centre):
        self.E1_centre = centre

    def set_E2_centre(self, centre):
        self.E2_centre = centre

    def set_CC_centre(self, centre):
        self.CC_centre = centre

    def set_data_field(self, field):
        self.data_field = gu.cut_field_for_fit(field)


    def run(self):
        self.current = 30.0*gu.get_polarity(self.data_field)

        self._setup_minuit_args()

        print 'Fitting to', len(self.data_field), 'data points'

        minimizer = minuit.Minuit(self._min_func, **self.minuit_kwargs)

        minimizer.migrad()

        return self._make_fit_dict(minimizer.values)


    def _setup_minuit_args(self):
        #This function is brutalistic but it makes the setup of minuit above look A LOT cleaner
        #setup CC dimension limits
        self.minuit_kwargs['CC_length'] = self.CC_length
        self.minuit_kwargs['limit_CC_length'] = (self.CC_length*0.95, self.CC_length*1.05)
        self.minuit_kwargs['CC_rInner'] = self.CC_rInner
        self.minuit_kwargs['limit_CC_rInner'] = (0.256, self.CC_rInner + 0.003)
        self.minuit_kwargs['CC_rOuter'] = self.CC_rOuter
        self.minuit_kwargs['limit_CC_rOuter'] = (self.CC_rOuter - 0.003, self.CC_rOuter + 0.003)
        #setup E1 dimension limits
        self.minuit_kwargs['E1_length'] = self.E1_length
        self.minuit_kwargs['limit_E1_length'] = (self.E1_length*0.95, self.E1_length*1.05)
        self.minuit_kwargs['E1_rInner'] = self.E1_rInner
        self.minuit_kwargs['limit_E1_rInner'] = (0.256, self.E1_rInner + 0.003)
        self.minuit_kwargs['E1_rOuter'] = self.E1_rOuter
        self.minuit_kwargs['limit_E1_rOuter'] = (self.E1_rOuter - 0.003, self.E1_rOuter + 0.003)
        #setup E2 dimension limits
        self.minuit_kwargs['E2_length'] = self.E2_length
        self.minuit_kwargs['limit_E2_length'] = (self.E2_length*0.95, self.E2_length*1.05)
        self.minuit_kwargs['E2_rInner'] = self.E2_rInner
        self.minuit_kwargs['limit_E2_rInner'] = (0.256, self.E2_rInner + 0.003)
        self.minuit_kwargs['E2_rOuter'] = self.E2_rOuter
        self.minuit_kwargs['limit_E2_rOuter'] = (self.E2_rOuter - 0.003, self.E2_rOuter + 0.003)
        #setup rest of parameters
        self.minuit_kwargs['CC_centre'] = self.CC_centre
        self.minuit_kwargs['E1_centre'] = self.E1_centre
        self.minuit_kwargs['E2_centre'] = self.E2_centre
        self.minuit_kwargs['thetaX'] = 0
        self.minuit_kwargs['limit_thetaX'] = (-0.01, 0.01)
        self.minuit_kwargs['thetaY'] = 0
        self.minuit_kwargs['limit_thetaY'] = (-0.01, 0.01)
        self.minuit_kwargs['px'] = 0
        self.minuit_kwargs['py'] = 0
        self.minuit_kwargs['pedantic'] = False


    def _min_func(self, CC_length, CC_rInner, CC_rOuter, E1_length, E1_rInner, E1_rOuter, \
                  E2_length, E2_rInner, E2_rOuter, CC_centre, E1_centre, E2_centre, thetaX, thetaY, \
                  px, py):

        E1_calc = mkfield.CalcField(E1_length, E1_rInner, E1_rOuter, self.E1_nLayers, self.E1_nTurns, \
                                    self.current, E1_centre, thetaX, thetaY, px, py)
        E2_calc = mkfield.CalcField(E2_length, E2_rInner, E2_rOuter, self.E2_nLayers, self.E2_nTurns, \
                                    self.current, E2_centre, thetaX, thetaY, px, py)
        CC_calc = mkfield.CalcField(CC_length, CC_rInner, CC_rOuter, self.CC_nLayers, self.CC_nTurns, \
                                    self.current, CC_centre, thetaX, thetaY, px, py)

        sum_br_res = 0
        sum_bp_res = 0
        sum_bz_res = 0
        
        for point in self.data_field:
            E1_Br, E1_Bp, E1_Bz = E1_calc.calc_field_at_point(point.r, point.phi, point.z)
            E2_Br, E2_Bp, E2_Bz = E2_calc.calc_field_at_point(point.r, point.phi, point.z)
            CC_Br, CC_Bp, CC_Bz = CC_calc.calc_field_at_point(point.r, point.phi, point.z)

            calc_Br = E1_Br + E2_Br + CC_Br
            calc_Bp = E1_Bp + E2_Bp + CC_Bp
            calc_Bz = E1_Bz + E2_Bz + CC_Bz
            
            sum_br_res += (point.Br - calc_Br)**2/(self.sigma*self.sigma)
            sum_bp_res += (point.Bphi - calc_Bp)**2/(self.sigma*self.sigma)
            sum_bz_res += (point.Bz - calc_Bz)**2/(self.sigma*self.sigma)

        chi_sq = sum_br_res + sum_bp_res + sum_bz_res
        print chi_sq/(len(self.data_field)-8)
        return chi_sq

        
    def _make_fit_dict(self, fitDict):
        coil_list = ['E1', 'CC', 'E2']

        ECE_fit_dict = {}

        for coil in coil_list:
            this_dict = {}
            #do easy ones first
            this_dict['px'] = fitDict['px']
            this_dict['py'] = fitDict['py']
            this_dict['thetaX'] = fitDict['thetaX']
            this_dict['thetaY'] = fitDict['thetaY']
            for key, value in fitDict.iteritems():
                if key.split('_')[0] == coil:
                    this_dict[key.split('_')[1]] = value
            this_dict['nTurns'] = self.nTurnsDict[coil]
            this_dict['nLayers'] = self.nLayersDict[coil]

            ECE_fit_dict[coil] = this_dict   

        return ECE_fit_dict
