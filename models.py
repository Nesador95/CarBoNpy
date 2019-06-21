# -*- coding: utf-8 -*-
"""
models.py - Temperature and Density Models

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

"""

import numpy as np

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
    T=T0 * (t / t0) ** (3 - 3 * gamma)

#    sup,sub=0,0
#    for name,index in idx.items():
#        if index!=99:
#            sup += mass[name]*y[index]*1.66054e-24
#            sub += y[index]
#    mu_gas = sup/sub
    ndens = dens * (t6 / t) ** 3   #dens/mu_gas*(t6/t)**3
    return T, ndens


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
    
    T = T0 * (t / t0) ** (3 - 3 * gamma)
    ndens = dens * (t6 / t) ** 3
    return T, ndens

def constantD(t,dens,T0 = 1.8e+4,gamma = 1.593, t0 = 100):
    '''
    This is a test of a model where the temperature decreases, but the number 
    density is constant.
    '''
    T=T0 * (t / t0) ** (3 - 3 * gamma)

    ndens = dens
    return T,ndens

def arrhenius(a,b,c,T,formula):
    '''
    Set the correct reaction rate formula for reactions contained in the 
    reactions file. More details for each rate formula are found at 
    http://kida.obs.u-bordeaux1.fr
    '''
    a *= 86400     # Convert time units from per second to per day

    # Formula choosing subroutine
    if formula==1:                    #Cosmic Ray Ionization
        zeta = 2.0e-17                #H2 ionization rate
        k = a * zeta

    elif formula==2:                  #Photodissociation (Draine)
        Av = 1                        #Av = visual extinction
        k = a * np.exp(-c * Av)

    elif formula==3:                  #Modified Arrhenius
        k = a * (T/300) ** b * np.exp(-c / T)

    elif formula==4:                  #ionpol 1   
        k = a * b * (0.62+0.4767 * c * np.sqrt(300. / T))

    elif formula==5:                  #ionpol 2
        k = a * b * (1 + 0.0967 * c * np.sqrt(300. / T) \
                     + (300 * c ** 2) / (10.526 * T))
    else:
        exit("Unknown formula, please check the input reactions file \
              for possible corruption")

    return k
