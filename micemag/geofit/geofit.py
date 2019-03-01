import sys
import os
import pickle

import scipy as sp
import numpy as np
import iminuit as minuit

import micemag.makefields.mkfieldclass as mkfield
import geofit_utils as gu


class GeoFit(object):
    def __init__(self, magParamDict, field=None, centre_pos=0, pc_change=10, current=30, sigma=0.002):
        self.pc_change = pc_change
        self.data_field = gu.cut_field_for_fit(field)
        self.centre_pos = centre_pos
        self.sigma = sigma
        
        if field != None:
            self.current = current*gu.get_polarity(field)
            if field[0].ID == 'Cartesian Data':
                print 'Need to pass field as polar coordinates to ', self.__class__.__name__
        else:
            self.current = 30    

        self.magParamDict = magParamDict
        self.setup_bracketing_fields(self.pc_change)
        

    def setup_bracketing_fields(self, pc_change):
        pc_mult_inc = 1.0 + (pc_change/100.0)
        pc_mult_dec = 1.0 - (pc_change/100.0)

        depth = self.magParamDict['rOuter'] - self.magParamDict['rInner']
        R_centre = self.magParamDict['rInner'] + depth/2.0
        #key: lt = longer/thinner (long in Z, thin in R)
        #     sf = shorter/fatter (short in Z, fat in R)
        lt_length = self.magParamDict['length']*pc_mult_inc
        lt_depth = depth*pc_mult_dec
        lt_rInner = R_centre - lt_depth/2.0
        lt_rOuter = R_centre + lt_depth/2.0

        self.lt_brfield = mkfield.CalcField(length=lt_length, rInner=lt_rInner, \
                                            rOuter=lt_rOuter, nLayers=self.magParamDict['nLayers'], \
                                            nTurns=self.magParamDict['nTurns'], current=self.current,\
                                            centrePos=self.centre_pos)

        sf_length = self.magParamDict['length']*pc_mult_dec
        sf_depth = depth*pc_mult_inc
        sf_rInner = R_centre - sf_depth/2.0
        sf_rOuter = R_centre + sf_depth/2.0

        self.sf_brfield = mkfield.CalcField(length=sf_length, rInner=sf_rInner, \
                                            rOuter=sf_rOuter, nLayers=self.magParamDict['nLayers'], \
                                            nTurns=self.magParamDict['nTurns'], current=self.current,\
                                            centrePos=self.centre_pos)

    def set_pc_change(self, pc_change):
        #This function is here incase anyone wants to play with how this affects stuff
        self.pc_change = pc_change

    def set_fit_data(self, field):
        self.data_field = gu.cut_field_for_fit(field)
        self.current = current*get_polarity(field)
    

    def run(self):
        if self.data_field == None:
            print 'No data field!'
            return

        minimizer = minuit.Minuit(self._min_func, mixing=0.5, limit_mixing=(0.35, 0.65), \
                                  scaling=1.0, limit_scaling=(0.9, 1.1), thetaX=0, \
                                  limit_thetaX=(-0.01, 0.01), thetaY=0, limit_thetaY=(-0.01, 0.01), \
                                  px=0, limit_px=(-0.1, 0.1), py=0, limit_py=(-0.1, 0.1), \
                                  pz=0, limit_pz=(-0.1, 0.1), pedantic=False)
        print 'Fitting to', len(self.data_field), 'data points'
        minimizer.migrad()

        return self._make_fit_dict(minimizer.values)

    def _min_func(self, mixing, scaling, thetaX, thetaY, px, py, pz):
        self.lt_brfield.set_angles_offsets(thetaX, thetaY, px, py)
        self.lt_brfield.set_centre_pos(self.centre_pos + pz)
        self.sf_brfield.set_angles_offsets(thetaX, thetaY, px, py)
        self.sf_brfield.set_centre_pos(self.centre_pos + pz)

        sum_br_res = 0
        sum_bp_res = 0
        sum_bz_res = 0

        
        
        for n, point in enumerate(self.data_field):
            sf_br, sf_bp, sf_bz = self.sf_brfield.calc_field_at_point(point.r, point.phi, point.z)
            lt_br, lt_bp, lt_bz = self.lt_brfield.calc_field_at_point(point.r, point.phi, point.z)
            #Which field gets multiplied by mixing and which gets multipied by mixing - 1 is important
            #to consider for when you want to reconstruct the fields from the parameters
            calc_Br = (sf_br*mixing + lt_br*(1-mixing))*scaling
            calc_Bp = (sf_bp*mixing + lt_bp*(1-mixing))*scaling
            calc_Bz = (sf_bz*mixing + lt_bz*(1-mixing))*scaling
            
            sum_br_res += (point.Br - calc_Br)**2/(self.sigma*self.sigma)
            sum_bp_res += (point.Bphi - calc_Bp)**2/(self.sigma*self.sigma)
            sum_bz_res += (point.Bz - calc_Bz)**2/(self.sigma*self.sigma)

        chi_sq = sum_br_res + sum_bp_res + sum_bz_res
        print chi_sq/(len(self.data_field)-7)
        return chi_sq

    def _make_fit_dict(self, fitDict):
        fitDict['pc_change'] = self.pc_change #not sure if i need anything else in here
        return fitDict

class GeoFitECE(object):
    def __init__(self):
        pass





class CalcGeoFit(object):
    def __init__(self):
        pass
