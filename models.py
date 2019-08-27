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
    t /= 86400

    T=T0 #* (t / t0) ** (3 - 3 * gamma)
    ndens = dens
    return T,ndens

def arrhenius(a,b,c,T,formula):
    '''
    Set the correct reaction rate formula for reactions contained in the 
    reactions file. More details for rate formulas 1-5 are found at 
    http://kida.obs.u-bordeaux1.fr
    
    Formula 6 is defined as f*K for grains. Van der Waals corrections are 
    added in Chemnet
    '''
    
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

def VdW(radius1,radius2,T,A,k=1.38e-23):
    """
    This function calculates the van der Waals correction factor for the coagulation 
    of grains of radius1 and radius2. 

    Input: Radii of grains 1 and 2, temperature, Hamaker constant. 
    Output: VdW correction coefficient. 

    We define the reduced Hamaker coefficient based on the work of Sarangi & 
    Cherchneff (2015) (who follow Sceats (1989)). See the derivation paperwork. 

    RT is the maximal real root of the polynomial equation from Sceats: 
    d/dr(V(r)-2kT ln(r))=0

    """
    a = radius1+radius2
    A_red = (1/(3*k*T))*A*radius1*radius2*a**4
    b = np.polynomial.Polynomial([-1*A_red,0,a**4,0,-2*a**2,0,1]).roots()
    RT = b[(b.imag==0) & (b.real>=0)].real.max()
    V = -(A/3)*radius1*radius2/a**2*(1/((RT/a)**2-1) + \
                                             1/(RT/a)**2 + \
                                             2*np.log(1-(a/RT)**2))
    W = (RT/a)**2*np.exp(-V/(k*T))

    return W
