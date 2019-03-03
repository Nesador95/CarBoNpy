# -*- coding: utf-8 -*-
"""
models.py - Temperature and Density Models

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

"""

##############################################################################
# The following is the Cherchneff Temperature and Density Model
##############################################################################

def cherchneffT(y, t, idx, mass, dens, T0 = 1.8e+4,
                gamma = 1.593, t0 = 100, t6 = 63.31):
    """
    Cherchneff & Dwek (2009) Temperature and Density Model
    quasi-adiabatic model based on Nozawa & Kozasa (2003)
    
    T(t) = T0*(t/t0)**(3(1-gamma))
    n(t) = rho(t)/mu_gas(t) 
    
    Default T0 at 100 days is 18,000K
    See C&D (2009) Table 1 for values for 20M_sun SNe
    mu_gas is the mean molecular weight of the gas.
    
    """ 
    
    t6 /= 365.25
    t0 /= 365.25  
    T=T0 * (t / t0) ** (3 - 3 * gamma)
#    sup,sub=0,0
#    for name,index in idx.items():
#        if index!=99:
#            sup += mass[name]*y[index]*1.66054e-24
#            sub += y[index]
#    mu_gas = sup/sub
    ndens = dens * (t6 / t) ** 3   #dens/mu_gas*(t6/t)**3
    return T, ndens

##############################################################################
# The following is the Yu Temperature and Density Model
###############################################################################

def YuT(t, dens, T0 = 3800, gamma = 4 / 3,
        t0 = 100, t6 = 63.31):
    """
    Yu et al. (2013) Temperature and Density Model
    "Basic Temperature Model"
    
    T(t) = 3.30e10 K.s /t 
    n(t) = <number>(t(6000k)/t)**3 cm**(-3) 
    
    Default T0 at 100 days is 3800K
    t(6000K) = 5.47e6 seconds = 63.31 days
    
    """
    
    t6 /= 365.25
    t0 /= 365.25
    T = T0 * (t / t0) ** (3 - 3 * gamma)
    ndens = dens * (t6 / t) ** 3
    return T, ndens
