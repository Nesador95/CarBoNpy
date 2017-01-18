# -*- coding: utf-8 -*-
'''
CarBon.py - Kinetic Molecular Chemistry Code for Supernova Ejecta

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

'''

from scipy.integrate import odeint
import datainput
import models
import numpy as np
from sys import exit

######
# Read data files
######

file_format,species_file,reactions_file,output_file,model_type,density = datainput.settings()

reactions,speciesidx,speciesmass,numspecies=datainput.kida_input(reac_file=reactions_file,spec_file=species_file)
count=len(reactions.index)

datainput.print_full(reactions)

#Functions that determine the reaction coefficients for each reaction. 
def arrhenius(a,b,c,T,formula):
    a *= 3.154e+7     # Convert units to year**-1

    # Formula choosing subroutine (from KIDA data set)
    if formula==1:                    #Cosmic Ray Ionization
        zeta = 2.0e-17                    #H2 ionization rate
        k = a*zeta

    elif formula==2:                  #Photodissociation (Draine)
        Av = 1                        #Av = visual extinction
        k = a*np.exp(-c*Av)

    elif formula==3:                  #Modified Arrhenius
        k = a*(T/300)**b*np.exp(-c/T)

    elif formula==4:                  #ionpol 1   
        k = a*b*(0.62+0.4767*c*np.sqrt(300./T))

    elif formula==5:                  #ionpol 2
        k = a*b*(1+0.0967*c*np.sqrt(300./T)+(300*c**2)/(10.526*T))
    else:
        exit("Unknown Formula!")

    return k

##Defining chemnet 

def chemnet(y,t):
    f=np.zeros(numspecies,float)
 
# Model Selection for T(t) and n(t)
    if model_type=='CD':
        T,Ndens=models.cherchneffT(y,t,speciesidx,speciesmass,Ndensinit)
    elif model_type=='Yu':
        T,Ndens=models.YuT(t,Ndensinit)
    else:
        exit("No Model Loaded, Exiting Now.")

    for num in range(count):        
        in1=reactions.loc[num,'Input1']
        in2=reactions.loc[num,'Input2']
        in3=reactions.loc[num,'Input3']
        out1=reactions.loc[num,'Output1']
        out2=reactions.loc[num,'Output2']
        out3=reactions.loc[num,'Output3']
        alpha=reactions.loc[num,'alpha']
        beta=reactions.loc[num,'beta']
        gamma=reactions.loc[num,'gamma']
        fo=reactions.loc[num,'Formula']
        if num==0:
            print('Still going! t={0}, Temp={1}'.format(t,T))

        if in2!=0 and in2!=99:
            f[in1] -= Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2]
            f[in2] -= Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2]
            f[out1] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2] 
            f[out2] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2]
            f[out3] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2]
        elif in2==99:
            f[in1] -= Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*(y[speciesidx['C']]+y[speciesidx['O']])
            f[out1] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*(y[speciesidx['C']]+y[speciesidx['O']]) 
            f[out2] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*(y[speciesidx['C']]+y[speciesidx['O']])
            f[out3] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*(y[speciesidx['C']]+y[speciesidx['O']])
        elif in2==0:
            f[in1] -= arrhenius(alpha,beta,gamma,T,fo)*y[in1]
            f[out1] += arrhenius(alpha,beta,gamma,T,fo)*y[in1] 
            f[out2] += arrhenius(alpha,beta,gamma,T,fo)*y[in1]
            f[out3] += arrhenius(alpha,beta,gamma,T,fo)*y[in1]
    return f



# Remember that time is measured in units of years.
time = np.linspace(60/365.25,5,1000000)

# Initial Values #########################################
yinit = np.zeros(numspecies,float)          
yinit[speciesidx['C']]  = 0.1            #Initial Carbon
yinit[speciesidx['O']]  = 1              #Initial Oxygen 
yinit[speciesidx['Si']] = 0            #Initial Silicon
Ndensinit = np.sum(yinit)*density
##########################################################

y = odeint(chemnet,yinit,time,mxstep=5000000,rtol=1e-13,atol=1e-13)

######
# Plot final abundances
######

abundance=[]

for i in range(1,numspecies):
    abundance.append(y[-1,i])

print(abundance)

######
# Write to data file
######

np.savez(output_file,time=time,speciesidx=speciesidx,y=y,abundance=abundance,speciesmass=speciesmass)
