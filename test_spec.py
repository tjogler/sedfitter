#!/usr/bin/env python
import os,sys,glob
import argparse
import utillities as ut
import SEDfitter as sed
import specShapes as spsh
from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
import testsed
import fitter_mod 
import constants as const
import crossSection as crs
import model as mod
import gamma_spectra2_mod as gs
import energyLoss as el

def plot(xlist,ylist,error_y=[],error_flaq=[0,0],marker_array=['a','b'],xtitle='E [MeV]',ytitle='EF',lpos='lower left'):
    '''
    plots the SED with individual colors for each experiment
    plots SED as points with errorbars (set error_flaq =1 for errors)
    plots models as solid curves (when error_flaq=0)
    adds legend with chi2 for models (set marker_array to item name)
    '''
    line=['-','--',':','-.']
    msymb=['o','s','D','v','^','h','*']
    fig=plt.figure()
    left,width=0.1,0.8
    rect1 = [left, 0.3, width, 0.6]
    rect2 = [left, 0.1, width, 0.2]
    ax1 = fig.add_axes(rect1)  #left, bottom, width, height
    ax2 = fig.add_axes(rect2, sharex=ax1)

    #plt.subplot(211)
    #ax=pyplot.subplots()
    counter=0
    counter_error=0
    counter_line=0
    
    for x in xlist:
        #plt.subplot(211)
        print 'adding plot %s of %s(%s)'%(counter,len(xlist),len(ylist))
        if np.size(x)!=np.size(ylist[counter]):
            print 'Error: x and y not the same dimensions (%i,%i) in point set %i'%(np.size(x),np.size(ylist[counter]),counter)
            exit()
        #print x,ylist[counter]
        if error_flaq[counter]==1:
            mstyle='o'
            #print 'error ' , error_y1
            #plt.loglog(x,ylist[counter],marker=msymb[counter_error],linestyle='None',label=marker_array[counter])
            ax1.errorbar(x,ylist[counter],marker=msymb[counter_error],yerr=error_y[counter_error],linestyle='None',label=marker_array[counter])
            #plt.yscale('log')
            #plt.xscale('log')
            counter_error+=1
        else:
            if counter_line >=4:
                counter_line=0
            ax1.plot(x,ylist[counter],linewidth=2.5,linestyle=line[counter_line],label=marker_array[counter])
            #plt.subplot(212)
            ax2.plot(x,abs(ylist[counter]-ylist[0])/ylist[0],linewidth=2.5,linestyle=line[counter_line],label=marker_array[counter])
            #plt.yscale('log')
            #plt.xscale('log')
            counter_line+=1
        counter+=1
    #calc min max of the y and x axis)

    counter=0
    counter_error=0
    for er in error_flaq:
        if er ==1:
            ylist.append(ylist[counter]+error_y[counter_error])
            ylist.append(ylist[counter]-error_y[counter_error])
            counter_error+=1
        counter+=1
    
    ymin=ut.min_listofarrays(ylist)*0.8
    ymax=ut.max_listofarrays(ylist)*1.2
    xmin=ut.min_listofarrays(xlist)*0.8
    xmax=ut.max_listofarrays(xlist)*1.2

    #plt.rc('text',usetex=True)
    #plt.subplot(211)
    if xmin >0:
        ax1.set_xscale('log')
        #(xscale='log')
        #plt.subplot(212)
        ax2.set_xscale('log')
    if ymin>0:
        #plt.subplot(211)
        ax1.set_yscale('log')

    print 'setting limits x=[%f,%f], y=[%f,%f]'%(xmin,xmax,ymin,ymax)
    ax1.set_ylim(ymin,ymax)
    ax1.set_ylabel(ytitle)
    ax2.set_xlabel(xtitle)
    ax1.set_xlim(xmin,xmax)
    ax1.legend(loc=lpos)
    #plt.subplot(212)
    ax2.set_ylabel(r'$\Delta\,$ '+ytitle)
    ax1.set_xlabel(xtitle)
    
    return fig



def plot_cross(typ=['ic','brems','pp']):
    Eg=np.logspace(6.,12.,300)
    Ee=np.logspace(7.,15.,100)
    Ne=spsh.plcut(Ee,[1.,9.,-3.,14.])
    Ne2=spsh.plcut(Ee,[1.,9.,-2.0,11.])
    ylist=[]
    xlist=[]
    taglist=[]
    elist=[]

    if 'brems' in typ:
        Ee/=10**9
        Eg/=10**9
        ylist.append(const.C*const.m2cm*crs.sigma_brems(Eg,Ee,Ne))
        print const.C*const.m2cm*crs.sigma_brems(Eg,Ee,Ne)
        #ylist.append(crs.sigma_brems(Eg,Ee,Ne2))
        ylist.append(const.C*const.m2cm*crs.sigma_brems_ee(Eg,Ee,Ne))
        #ylist.append(crs.sigma_brems_ee(Eg,Ee,Ne2))
        ylist.append(const.C*const.m2cm*gs.sigmaB_test(Ne,Eg,Ee))
        xlist.append(Eg)
        xlist.append(Eg)
        xlist.append(Eg)
        #xlist.append(Eg)
        #xlist.append(Eg)
        for t in['brems_N -1.5','brems_ee -1.5','brems dima']: 
            taglist.append(t)
            elist.append(0)
    if 'ic' in typ:
        Ee/=10**9
        Eg/=10**9
        #print Ee,Eg
        source=sed.source()
        source.set_radiation_densities(0.25,0.84,3.0e-3,0.9)#W51C densities Abdo et al 2009
        source.set_radiation_photons(-8.,4.,0.1)
        '''
        print 'CMB: ',source.phCMB
        print 'Dust : ',source.phDust
        print 'Star: ',source.phStar
        print 'tot: ',source.phTot
        '''
        for ph in [source.phTot,source.phCMB,source.phDust,source.phStar]:
            def sigma_fcn(x):
                sigma=crs.sigma_ic(x,source.Eph,Ee)
                #print 'sigma: ',np.size(sigma)
                return const.C*const.m2cm*np.dot(ph,np.dot(sigma,np.ones_like(Ee)))

            sigma=np.frompyfunc(sigma_fcn,1,1)
            ylist.append(sigma(Eg))
            xlist.append(Eg)
            elist.append(0)
        for t in['IC total','IC CMB','IC Dust','IC Star']: 
            taglist.append(t)
        plot([source.Eph,source.Eph,source.Eph,source.Eph],[source.phTot,source.phCMB,source.phDust,source.phStar],[],elist,['Tot','CMB','Dust','Star'],r'E [eV]',r'$N_{ph}$')
    
    #print ylist
    
    plot(xlist,ylist,[],elist,taglist,r'E [eV]',r'$N_{e}\times \sigma_{\mathrm{brems}}$')

    #plot([source.Eph,source.Eph,source.Eph,source.Eph],[source.phTot,source.phCMB,source.phDust,source.phStar],[],elist,['Tot','CMB','Dust','Star'],r'E [eV]',r'$N_{ph}$')
    
    plt.show()



def plot_spectrum(typ=['ic','brems','pp']):
    '''plots spectra'''
    Eg=np.logspace(6.,12.,100)
    Ee=np.logspace(6.,15.,1000)
    parameter=[1.,9.,-2.5,10.]
    #Ne2=spsh.bplcut(Ee,[49.4,-1.5,-2.9,10.176,13.])
    Ne2=spsh.sbplcut(Ee,[49.4,9.,-1.5,10.301,2.3,13.])
    #print Ne
    Ne=spsh.sbplcuty(Ee,[49.4,9.,-1.5,10.301,2.3,13.])#W51C values abdo et al 2009
    ENe=Ee*Ne2
    spec=sed.spectrum(spsh.sbplcut,[49.4,9.,-1.5,10.301,2.3,13.],1.e6,1.e15,1000.)
    
    #Ne2=spsh.plcut()
    ylist=[]
    xlist=[]
    taglist=[]
    elist=[]
    source=sed.source()
    source.set_radiation_densities(0.25,0.84,3.0e-3,0.9)#W51C densities Abdo et al 2009
    source.set_radiation_photons(-8.,4.,0.1)
    source.set_am_prop(10.,10.,1.)
    ppcrosspath='/Users/jogler/Physik/Fermi/scripts/pp_data/'
    ppcross=crs.crossSection(ppcrosspath)

    if 'brems' in typ:
        brems=mod.model(spec,'brems',[],source,ppcross)
        En=brems.spec.get_xvals()*brems.spec.value()
        print En-ENe
        ylist.append((Eg/1.e9)*brems.Brems_spec(ENe,Ee/1.e9,Eg/1.e9))
        #ylist.append(brems.Brems_spec(ENe,Ee/10**9,Eg/10**9))
        # ylist.append(brems.Brems_spec_ee(Ne,Ee,Eg))
        xlist.append(Eg)
        #xlist.append(Eg)
        taglist.append('Brems N')
        #taglist.append('Brems ee')
        elist.append(0)
        #elist.append(0)

    if 'ic' in typ:
        for ph in [source.phTot,source.phCMB,source.phDust,source.phStar]:
            ic=mod.model(Ne,'ic',[],source,ppcrosspath)
            ylist.append(Eg/1.e9*ic.IC_spec_gamma(ph/source.Eph,source.Eph,ENe,Ee/10**9,Eg/10**9))
            xlist.append(Eg)
            elist.append(0)
            
        for t in['IC total','IC CMB','IC Dust','IC Star']: 
            taglist.append(t)
    
    if 'pp' in typ:
        pp=mod.model(Ne,'pp',[],source,ppcross)
        ylist.append(pp.EdQdE_pp(Ne,Ee/10**9,Eg/10**9))
        xlist.append(Eg)
        elist.append(0)
        taglist.append(r'$\mathrm{pp}_{\gamma} $')
        ''' plot e+ e- sec spec'''
        ylist.append(pp.EdQdE_pp(Ne,Ee/10**9,Eg/10**9,1))
        xlist.append(Eg)
        elist.append(0)
        taglist.append('pp_e-')
        ylist.append(pp.EdQdE_pp(Ne,Ee/10**9,Eg/10**9,2))
        xlist.append(Eg)
        elist.append(0)
        taglist.append('pp_e+')
    if 'comp' in typ:
        if 'brems' in typ:
            print 'nH: %.2f nHeff: %.2f Zeff: %.2f'%(source.nH,source.nHeff,source.Zeff)
            bremsSp=gs.brems_spectrum(ENe,Ee/10**9,source.nHeff)
            ylist.append(Eg/10**9*bremsSp(Eg/10**9))
            xlist.append(Eg)
            elist.append(0)
            taglist.append('brems dima')
        if 'ic' in typ:
            icSp=gs.IC_spectrum(source.phTot/source.Eph,source.Eph,ENe,Ee/10**9)
            ylist.append(Eg/10**9*icSp(Eg/10**9))
            xlist.append(Eg)
            elist.append(0)
            taglist.append('ic dima')
        if 'pp' in typ:
            ppSp=gs.EdQdE_pp(Ne,Ee/10**9,source.nH,0)
            ylist.append(ppSp(Eg/10**9))
            xlist.append(Eg)
            elist.append(0)
            taglist.append('pp dima')
        
    fig0=plot(xlist,ylist,[],elist,taglist,r'$\mathrm{E}_{\gamma}\, \mathrm{ [eV]}$',r'$E^2F_{\gamma}\, [\mathrm{erg}\,\mathrm{ cm}^{-3} \mathrm{s}^{-1}]$')
    #plt.show(fig0)

    xlist=[]
    ylist=[]
    elist=[]
    taglist=[]
    
    xlist.append(Ee)
    xlist.append(Ee)
    ylist.append(Ne)
    ylist.append(Ne2)
    elist.append(0)
    elist.append(0)
    taglist.append('Ne')
    taglist.append('Ne2')
    fig1=plot(xlist,ylist,[],elist,taglist,r'$\mathrm{E}_{\gamma}\, \mathrm{ [eV]}$',r'$N_e\, [\mathrm{N}_{e}\mathrm{ cm}^{-3} \mathrm{s}^{-1}]$')
    plt.show()

def plot_sync_loss( B=np.logspace(0.,2.,6),Ee=np.logspace(7.,12.,100)/1.e9,pitch=np.linspace(0.1,1.,5)):
    
    #Ee=np.logspace(7.,12.,100)/1.e9
    #B=np.array(1.,10.,100.)
    #pitch=np.linspace(0.1,1.,5)
    #B=np.logspace(0.,2.,6)
    print B
    xlist=[]
    ylist=[]
    ylist2=[]
    ylist3=[]
    xlist3=[]
    elist3=[]
    taglist3=[]
    elist=[]
    taglist=[]

    spec=el.energyLoss()

    #print spec.Sync_Edot(B,Ee,pitch[4])

    for p in pitch:
        #print B,Ee,p
        ylist.append(spec.Sync_Edot(B[0],Ee,p))
        ylist2.append(Ee/spec.Sync_Edot(B[0],Ee,p)/const.yr2s)
        xlist.append(Ee)
        elist.append(0)
        taglist.append(r'$\sin( \alpha)$=%.2f'%p)
    
    for b in B:
        ylist3.append(Ee/spec.Sync_Edot(b,Ee,1.)/const.yr2s)
        taglist3.append(r'$B=%i \mu\mathrm{G}$'%b)
        xlist3.append(Ee)
        elist3.append(0)

    plot(xlist,ylist,[],elist,taglist,r'$\mathrm{E}\, \mathrm{ [GeV]}$',r'$\mathrm{dE}/\mathrm{d}t\, \mathrm{ [GeV/s]}$','lower right')
    plt.show()
    plot(xlist,ylist2,[],elist,taglist,r'$\mathrm{E}\, \mathrm{ [GeV]}$','cooling time  [yr]','lower left')
    plt.show()
    plot(xlist3,ylist3,[],elist3,taglist3,r'$\mathrm{E}\, \mathrm{ [GeV]}$','cooling time  [yr]','lower left')
    plt.show()



def plot_sec_spec():
    particlePoints=1000
    Eg=np.logspace(7.,13.,100)/1.e9
    Ee=np.logspace(6.,15.,particlePoints)/1.e9
    parameter=[1.,9.,-2.5,10.]
    #Ne=spsh.bplcut(Ee,[49.4,-2.43,-2.7,1.69,5.])
    Ne=spsh.sbplcuty(Ee,[49.4,1.,-1.5,1.176,-1.4,5.])
    
    #Ne2=spsh.plcut()
    ylist=[]
    xlist=[]
    taglist=[]
    elist=[]
    ppcrosspath='/Users/jogler/Physik/Fermi/scripts/pp_data/'
    ppcross=crs.crossSection(ppcrosspath)

    source=sed.source(6000,30.,3.e4)
    source.set_radiation_densities(0.25,0.84,3.0e-3,0.9)#W51C densities Abdo et al 2009
    source.set_radiation_photons(-8.,4.,0.1)
    source.set_am_prop(10.,10.,1.)
    
    fluxScale=Eg/(4.*np.pi*(source.dist*const.PARSEC*const.m2cm)**2.)/const.erg2GeV
    print 'Scale: ',fluxScale
    pp=mod.model(Ne,'pp',[],source,ppcross)
    ppSpec=pp.nucFactor*pp.EdQdE_pp(Ne,Ee,Eg)
    #secEPoints=np.logspace(np.log10(Eg[0]),np.log10(Eg[len(Eg)-1]),particlePoints)
    eSecSpec=pp.nucFactor*(pp.EdQdE_pp(Ne,Ee,Ee,1)+pp.EdQdE_pp(Ne,Ee,Ee,2))*source.age*const.yr2s
    print 'points in secondary spec: ',np.size(eSecSpec)
    
    secBrems=pp.Brems_spec(eSecSpec,Ee,Eg)
    secIC=pp.IC_spec_gamma(source.phTot,source.Eph,eSecSpec,Ee,Eg)

    ylist.append(ppSpec*fluxScale)
    ylist.append(secBrems*fluxScale)
    ylist.append(secIC*fluxScale)
    
    xlist.append(Eg)
    xlist.append(Eg)
    xlist.append(Eg)
    elist.append(0)
    elist.append(0)
    elist.append(0)
    taglist.append('pp')
    taglist.append('brems sec')
    taglist.append('IC sec')
    
    fig=plot(xlist,ylist,[],elist,taglist,r'$\mathrm{E}_{\gamma}\, \mathrm{ [GeV]}$',r'$EF_{\gamma}\, [\mathrm{erg}\, \mathrm{ cm}^{-2} \mathrm{s}^{-1}]$')
    plt.show(fig)
