# -*- coding: utf-8 -*-
"""
models.py - Temperature and Density Models

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

"""

"""
Cherchneff & Dwek (2009) Temperature and Density Model
quasi-adiabatic model based on Nozawa & Kozasa (2003)

T(t) = T0*(t/t0)**(3(1-gamma))
n(t) = rho(t)/mu_gas(t) 

Default T0 at 100 days is 18,000K
See C&D (2009) Table 1 for values for 20M_sun SNe
mu_gas is the mean molecular weight of the gas.
""" 

def cherchneffT(y,t,idx,mass,dens,T0=1.8e+4,gamma=1.593,t0=100,t6=63.31):
    t6 /= 365.25
    t0 /= 365.25  
    T=T0*(t/t0)**(3-3*gamma)
    sup,sub=0,0
    for name,index in idx.items():
        if index!=99:
            sup += mass[name]*y[index]
            sub += y[index]
    mu_gas = sup/sub
    ndens = dens/mu_gas*(t0/t)**3
    return T,ndens

"""
Yu et al. (2013) Temperature and Density Model
"Basic Temperature Model"

T(t) = 3.30e10 K.s /t 
n(t) = <number>(t(6000k)/t)**3 cm**(-3) 

Default T0 at 100 days is 3800K
t(6000K) = 5.47e6 seconds = 63.31 days
"""

def YuT(t,dens,T0=3800,t0=100,t6=63.31):
    t6 /= 365.25
    t0 /= 365.25
    T = T0*(t0/t)
    ndens = dens*(t6/t)**3
    return T,ndens
