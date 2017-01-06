#!/usr/bin/python
#title			:CarBon.py
#description		:Kinetic Molecular Chemistry Code for Supernova Ejecta
#author			:Ethan Deneault
#date			:20161027
#version		:3.6.0
#usage			:python CarBon.py
#notes			:See UPDATES
#python_version		:3.5.2
#==============================================================================
# 
#IMPORTANT:
#--All KIDA Data Sets Must have a header line included in it, may need to put it in / edit yourself! 
#--If there are 3 output species, make sure that at least ONE is placed on the first data line. 
#Please Read data/kida_readme_2016-x-x_1.dat for info

from scipy.integrate import odeint
import pandas as pd
import numpy as np
from sys import exit

######
# Input Files
######

species_file='data/kida_spec_C_O_Si_only.dat'
reactions_file='data/kida_reac_C_O_Si_only.dat'

######
# Output Data File
######

output_file='output/D_2017/CD_Model/output_C1_O0.01_Si1.dat'

#Read species file, create 2 dictionaries 
#speciesidx relating species name to index 
#speciesmass relating species name to 'mass' 

speciesidx={}
speciesmass={}
with open(species_file,'r') as my_file:
    for line in my_file:
        columns = line.strip().split()
        speciesidx[columns[0]]=int(columns[24])
        sum=0
        for i in range(2,24):
            num = int(columns[i])
            sum += num
        speciesmass[columns[0]]=sum 

numspecies=len(speciesidx)+1 # Total number of chemical species

#Convert Inputs to Species Number
#Some inputs are missing (typically Output2 or Output3) Don't change those values. Let Pandas do that. 

convert = lambda x: int(speciesidx[x]) if x!='' else None

#Read the reactions file and replace species names with index numbers from species file.
#First two lines read the file with properly formatted columns. See: http://stackoverflow.com/questions/40663510/pandas-read-fwf-ignoring-values
#Make sure that headers are set correctly to cover negative signs!
#Third line replaces NaN with species number zero
#Fourth line forces index number columns to be integers. 

reactions = pd.read_fwf(reactions_file, header=None, skiprows=1,comment='#',converters={0:convert,1:convert,2:convert,3:convert,4:convert})
reactions.columns = pd.read_csv(reactions_file, delim_whitespace=True,nrows=1).columns
reactions.fillna(0,inplace=True)
reactions[['Input1','Input2','Output1','Output2','Output3','Re']] = reactions[['Input1','Input2','Output1','Output2','Output3','Re']].astype(int)

count=len(reactions.index)

#### Test To see if Reactions File is read correctly
def print_full(x):
    pd.set_option('display.max_rows', len(x))
    pd.set_option('display.max_columns', 1000)
    pd.set_option('display.expand_frame_repr', False)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.expand_frame_repr')

print_full(reactions)
#for name,index in speciesidx.items():
#    if index != 0:
#        print(name,index)
#exit()
####

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

stepnum = 0

def chemnet(y,t):
    f=np.zeros(numspecies,float)
    # Define Temperature (All data from Cherchneff and Dwek 2009)  
    t0 = 100/365.25   # 0.273... years is approximately 100 days
    T0 = 1.8e+4       # Temperature at 100 days (1.8e+4 from C&D, 6000 from Yu et al)
    gamma = 1.593     # "adiabatic" index
    T=T0*(t/t0)**(3-3*gamma)
    # "Basic temperature model" From Yu et al:
    #T0 = 3.3e+10/3.154e+07
    #t0 = 63.6929/365.25
    #T = T0/t
    
    ##Test Code for number density (see Deneault et al 06 and Yu et al 13)
    #Ndensinit = 1.1e+10
    Ndens = Ndensinit*(t0/t)**3 

    for num in range(count):        
        in1=reactions.loc[num,'Input1']
        in2=reactions.loc[num,'Input2']
        out1=reactions.loc[num,'Output1']
        out2=reactions.loc[num,'Output2']
        out3=reactions.loc[num,'Output3']
        alpha=reactions.loc[num,'alpha']
        beta=reactions.loc[num,'beta']
        gamma=reactions.loc[num,'gamma']
        fo=reactions.loc[num,'Fo']
        if num==0:
            print('Still going! t={0}, Temp={1}'.format(t,T))

        #print(in1,in2,out1,out2,out3,alpha,beta,gamma,fo)
        #print(arrhenius(alpha,beta,gamma,t,fo)) 
        if in2 != 0:
            f[in1] -= Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2]
            f[in2] -= Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2]
            f[out1] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2] 
            f[out2] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2]
            f[out3] += Ndens*arrhenius(alpha,beta,gamma,T,fo)*y[in1]*y[in2]
        elif in2 == 0:
            f[in1] -= arrhenius(alpha,beta,gamma,T,fo)*y[in1]
            f[out1] += arrhenius(alpha,beta,gamma,T,fo)*y[in1] 
            f[out2] += arrhenius(alpha,beta,gamma,T,fo)*y[in1]
            f[out3] += arrhenius(alpha,beta,gamma,T,fo)*y[in1]
    return f



# Remember that time is measured in units of years.
time = np.linspace(60/365.25,5,1000000)

# Initial Values #########################################
yinit = np.zeros(numspecies,float)          
yinit[speciesidx['C']]  = 1            #Initial Carbon
yinit[speciesidx['O']]  = 0.01              #Initial Oxygen 
yinit[speciesidx['Si']] = 1            #Initial Silicon
Ndensinit = np.sum(yinit)*1e10
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

np.savez(output_file,time=time,speciesidx=speciesidx,y=y,abundance=abundance)
