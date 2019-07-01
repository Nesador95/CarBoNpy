# -*- coding: utf-8 -*-
'''
CarBon_1.0.1.py - Kinetic Molecular Chemistry Code for Supernova Ejecta

Copyright (c) 2012-2016 Ethan Deneault 
This file is part of CarBoN

''' 
import numpy as np
import matplotlib.pyplot as plt
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode
from sys import exit


import datainput as d
import models as m


"""
Import data from the settings file, the KIDA database files and the initial abundances file. 
"""

file_format, species_file, reactions_file, output_file, model_type, density, Tinit, start_time, end_time, outfile = d.settings()

kida_reac,kida_spec,spec_dict = d.KIDA_input(reactions_file,species_file)

abund_df = d.abundances(spec_dict)


"""
This list contains the monoatomic molecules that result from the use of 'M' 
in the KIDA files
C, O, Si
"""

monoatomic_list = [1, 2, 22]


count=len(list(kida_reac.index))



def chemnet(t,y):
    '''
    This function defines the rhs of the system of differential equations used in the chemical 
    network. It dynamically builds the equations from reading the reactions file, and defines the Arrhenius 
    coefficients from the coefficients given per reaction. The temperature model_type must be assigned to 
    either 'CD' or 'Yu'. 
    '''

    if model_type=='CD':
        T,Ndens=m.cherchneffT(y,t,kida_spec['species_num'].to_dict(),
                                   kida_spec['atom_num'].to_dict(),
                                   dens=Ndensinit,
                                   T0=Tinit)
    elif model_type=='Yu':
        T,Ndens=m.YuT(t,dens=Ndensinit,T0=Tinit)
    elif model_type=='Cons':
        T,Ndens=m.constantD(t,dens=Ndensinit,T0=Tinit)
    else:
        exit("No Model Loaded, Exiting Now.")

    f=np.zeros([len(kida_spec.index)]) # Define the rhs array

    #f -= 3*y/t


    for num in range(count):
        
        in1=kida_reac.loc[num]['Input1']
        in2=kida_reac.loc[num]['Input2']
        out1=kida_reac.loc[num]['Output1']
        out2=kida_reac.loc[num]['Output2']
        out3=kida_reac.loc[num]['Output3']     
        alpha=kida_reac.loc[num]['alpha']
        beta=kida_reac.loc[num]['beta']
        gamma=kida_reac.loc[num]['gamma']
        fo=kida_reac.loc[num]['Fo']

        if num==0:
            print('Still going! t={0}, Temp={1}'.format(t,T))

        if in2!=0 and in2!=99:
            f[in1]  -=  m.arrhenius(alpha,beta,gamma,T,fo) * y[in1] * y[in2]
            f[in2]  -=  m.arrhenius(alpha,beta,gamma,T,fo) * y[in1] * y[in2]
            f[out1] +=  m.arrhenius(alpha,beta,gamma,T,fo) * y[in1] * y[in2]           
            if np.isnan(out2) == False and out2!=0:
                f[out2] += m.arrhenius(alpha,beta,gamma,T,fo) * y[in1] * y[in2]
            if np.isnan(out3) == False and out3!=0:
                f[out3] += m.arrhenius(alpha,beta,gamma,T,fo) * y[in1] * y[in2]  

        elif in2==0:
            f[in1] -= m.arrhenius(alpha, beta, gamma, T, fo) * y[in1]
            f[out1] += m.arrhenius(alpha, beta, gamma, T, fo) * y[in1]
            if np.isnan(out2) == False and out2!=0:
                f[out2] += m.arrhenius(alpha, beta, gamma, T, fo) * y[in1]
            if np.isnan(out3) == False and out3!=0:
                f[out3] += m.arrhenius(alpha, beta, gamma, T, fo) * y[in1]

########### Reactions with moderators needs rewrite!
# 
#        elif in2==99:
#            f[in1] -=  Ndens*m.arrhenius(alpha,beta,gamma,T,fo) \
#                      * y[in1] * (kida_spec.loc[0]['species_num']\
#                         + kida_spec.loc[2]['species_num'] \
#                      + monoatomic_list[2])
#            f[out1] +=  Ndens*m.arrhenius(alpha,beta,gamma,T,fo) \
#                       * y[in1] * (kida_spec.loc[0]['species_num'] + kida_spec.loc[1]['species_#'] \
#                       + monoatomic_list[2])
#            if isinstance(out2,int) == True:
#                f[out2] += Ndens*m.arrhenius(alpha,beta,gamma,T,fo) * y[in1] \
#                       * (kida_spec.loc[0]['species_num'] + kida_spec.loc[1]['species_#'] \
#                       + monoatomic_list[2])
#            if isinstance(out3,int) == True:
#                f[out3] += Ndens*m.arrhenius(alpha,beta,gamma,T,fo) * y[in1] \
#                        * (kida_spec.loc[0]['species_num'] + kida_spec.loc[1]['species_#'] \
#                       + monoatomic_list[2])
############ NEEDS CLEANUP

    return f

'''
Initialize the variables
'''

yinit = np.zeros([len(kida_spec.index)])     
yinit[abund_df.Species.values] = abund_df.Abundance.values 

print(yinit)

Ndensinit = density #np.sum(yinit)*density

#testchem=chemnet(start_time,yinit)

#exit()

'''
Initialize Assimulo, calculate with CVode.
'''

start_time *= 86400
end_time *= 86400

model=Explicit_Problem(chemnet,yinit,start_time)
model.name='Chemnet Test'

sim=CVode(model)

sim.atol=1.e-12
sim.rtol=1.e-12
sim.maxord=3
sim.discr='BDF'
sim.iter='Newton'

t,y=sim.simulate(end_time)

#sim.plot()

plt.ylim([1e-2,2e8])
plt.semilogy(t,y[:,1],label='C')
plt.semilogy(t,y[:,2],label='O')
plt.semilogy(t,y[:,3],label='C2')
plt.semilogy(t,y[:,4],label='CO')
plt.semilogy(t,y[:,5],label='C3')
plt.semilogy(t,y[:,6],label='C4')
plt.semilogy(t,y[:,7],label='C5')
plt.semilogy(t,y[:,8],label='C6')
plt.semilogy(t,y[:,9],label='C7')
plt.semilogy(t,y[:,10],label='C8')
plt.semilogy(t,y[:,11],label='C9')
plt.semilogy(t,y[:,12],label='C10')

plt.legend()
plt.show()

##############################################################################
# The following writes to a data file
##############################################################################
#
#np.savez(output_file,time=time,speciesidx=kida_spec['species_num'].to_dict(), y = y,
#         abundance=abundance,speciesmass=kida_spec['atoms_num'].to_dict())
