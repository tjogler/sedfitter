#!/usr/bin/env python
import scipy as sp
import numpy as np
import constants as const
import fitter
import specShapes as spsh
import crossSection as crss
import SEDfitteer as sed
import random

class testdata:
    def __init__(self,npoints,ipara,emin,emax,specshape,errorsize):
        self.npts=npoints
        self.xmin=emin
        self.xmax=emax
        self.shape=specshape
        self.error=errorsize
        self.para=ipara
        _make_data

    def _make_data(self):
        print 'Creating test data of %i points using shape %s and parameters %s'%(self.npts,self.shape,self.para)
        print 'points are modulated with gaussian distribution with sigma=%f'%self.error
        self.spectrum=sed.spectrum(spsh.name(self.shape),self.para,self.xmin,self.xmax,self.npts)
        print 'test ',self.spectrum.par
        yvals=self.spectrum.value()
        '''make gaussian errorbars around the yvals that are gaussian randomized where self.error= sigma '''
        if random.uniform(0,1)>0.5:
            sig=1
        else:
            sig=-1
        
        output=yvals+(sig*np.random.normal(yvals,self.error*np.ones_like(yvals)))
        return output
        
