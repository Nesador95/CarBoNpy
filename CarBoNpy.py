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


def chemnet(t,y):
    '''
    This function defines the rhs of the system of differential equations 
    used in the chemical network. It dynamically builds the equations from 
    the DataFrame constructed by m.KIDA_Input()
    '''

    if model_type=='Cons':
        T,Ndens=m.constantD(t,dens=Ndensinit,T0=Tinit)
    else:
        exit("No Model Loaded, Exiting Now.")

    f=np.zeros([len(kida_spec.index)+len(grain_spec.index)]) # Define the rhs array

    f -= 3*y/t


    for num in range(len(list(kida_reac.index))):
        
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

        
        elif in2!=0 and in2!=99:
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

    for num in range(len(list(grains_reac.index))):

        in1=grains_reac.loc[num]['Input1']
        in2=grains_reac.loc[num]['Input2']
        out1=grains_reac.loc[num]['Output1']
        out2=grains_reac.loc[num]['Output2']     
        fijk=kida_reac.loc[num]['f_ijk']
        kij=kida_reac.loc[num]['K_ij']
        Hamaker=kida_reac.loc[num]['Hamaker']
        fo=kida_reac.loc[num]['Fo']

        f[in1] -= kij * np.sqrt(T) * m.VdW(r1,r2,T,Hamaker)
        f[in2] -= kij * np.sqrt(T) * m.VdW(r1,r2,T,Hamaker)
        f[out1] += fijk * kij * np.sqrt(T) * m.VdW(r1,r2,T,Hamaker)
        f[out2] += (1-fijk) * kij * np.sqrt(T) * m.VdW(r1,r2,T,Hamaker)



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


"""
Import data from the settings file, the KIDA database files and the initial abundances file. 
"""

# This is to read files from whatever they are stored in. 


file_format, MODEL_FILE, model_type, density, Tinit, start_time, end_time, outfile = d.settings()

kida_reac, kida_spec, spec_dict, grains_reac = READ?(MODEL_FILE)

abund_df = d.abundances(spec_dict)

"""
Initialize the variables
"""

yinit = np.zeros([len(kida_spec.index)])     
yinit[abund_df.Species.values] = abund_df.Abundance.values 

print(yinit)

Ndensinit = density #np.sum(yinit)*density

"""
Initialize Assimulo, calculate with CVode.
"""


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

t = np.array(t)/86400
plt.ylim([1e1,2e8])
plt.semilogy(t,y[:,1],label='C')
plt.semilogy(t,y[:,2],label='O')
plt.semilogy(t,y[:,3],label='C2')
plt.semilogy(t,y[:,4],label='CO')
plt.semilogy(t,y[:,5],label='O2')
plt.semilogy(t,y[:,6],label='C3')
plt.semilogy(t,y[:,7],label='C4')


plt.legend()
plt.show()

##############################################################################
# The following writes to a data file
##############################################################################
#
#np.savez(output_file,time=time,speciesidx=kida_spec['species_num'].to_dict(), y = y,
#         abundance=abundance,speciesmass=kida_spec['atoms_num'].to_dict())
