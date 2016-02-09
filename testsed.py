#!/usr/bin/env python
import scipy as sp
import numpy as np
import constants as const
import fitter
import specShapes as spsh
import crossSection as crss
import SEDfitter as sed
import random

class testdata:
    def __init__(self,npoints,ipara,emin,emax,specshape,errorsize):
        self.npts=npoints
        self.xmin=emin
        self.xmax=emax
        self.shape=specshape
        self.error=errorsize
        self.para=ipara
        

    def make_data(self):
        print 'Creating test data of %i points using shape %s and parameters %s'%(self.npts,self.shape,self.para)
        print 'points are modulated with gaussian distribution with sigma=%f'%self.error
        
        spectrum=sed.spectrum(spsh.func(self.shape),self.para,self.xmin,self.xmax,self.npts)
        '''Make the y data from the spectrum'''
        yvals=spectrum.value()

        '''make gaussian errorbars around the yvals that are gaussian randomized where self.error= sigma '''
        output=np.random.normal(yvals,self.error*yvals)
        xout=np.logspace(np.log10(self.xmin),np.log10(self.xmax),self.npts)
        #print xout
        return np.array([xout,output,np.zeros_like(xout),self.error*yvals,self.error*yvals,5*np.ones_like(xout)])
