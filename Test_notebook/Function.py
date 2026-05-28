#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 22:40:56 2024

@author: mokka
"""
import numpy as np
import random


def modelG(z, m_halo, cospar, params):

    '''
    DM Model

    Input:
    ----------
        z(redshift): float,
        m_halo: float, halo mass in M_sun/h
        cospar: list with the cosmological parameters
        params: list with the BFC parameters
        
    Output
    ----------
        object host
    '''

    import baryonification as bfc
    from scipy.interpolate import splrep,splev
    import pyccl as ccl

    par = bfc.par()

    par.cosmo.Om = cospar[0]
    par.cosmo.Ob = cospar[1]
    par.cosmo.s8 = cospar[2]
    par.cosmo.ns = cospar[3]
    par.cosmo.h0 = cospar[4]
    par.cosmo.z  = z

    par.baryon.eta = params[0]  #0.1  # slope total stellar to halo fraction
    par.baryon.deta = params[1] #0.22  # slope central stellar to halo fraction
    par.baryon.Nstar = params[2]#0.028  # Normalisation parameter
    par.baryon.ciga = params[3] #0.1 # Normalisation inner cold gas fraction (not important for LSS)
    par.baryon.thco = params[4] #0.3 # Core radius
    par.baryon.Mc = params[5]   #10**11.8609661 # Halo mass scale at which the slipe of hga is -1.5
    par.baryon.mu = params[6]   #1.04918393 # slope rate variation with halo mass
    par.baryon.delta = params[7]#5.27091664 # slope of the outer truncation of hga

    # Parameters which are always fixed
    par.baryon.gamma = params[8]
    par.baryon.alpha = 1.0
    par.baryon.a_nth = 0.001
    par.code.beta_model = 1.0
    par.code.eps0 = 4.0
    par.code.eps1 = 0.5
    par.code.AC_model = 5.0
    par.code.q0 = 0.075
    par.code.q1 = 0.25
    par.code.q2 = 0.7

    #m_halo =  #in M_sun/h

    rbin = np.geomspace(1e-4, 500, 100) #in Mpc/h
    host = bfc.frb_host_pdf(par,rbin,m_halo,rvir_mult=5)
    return host
