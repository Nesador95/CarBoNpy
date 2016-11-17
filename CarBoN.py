#!/usr/bin/python
#title			:CarBon.py
#description		:Kinetic Molecular Chemistry Code for Supernova Ejecta
#author			:Ethan Deneault
#date			:20161027
#version		:3.1.5
#usage			:python CarBon.py
#notes			:See UPDATES
#python_version		:3.5.1
#==============================================================================
#assoc.dat contains all reactions: A+B->C with columns set as:  
#[input 1,input 2,output,alpha,beta,gamma] 
#
#dblassoc.dat contains all reactions: A+B->C+D with columns set as:
#[input 1,input 2,output 1,output 2,alpha,beta,gamma]
#
#breakup.dat contains all reactions: A->B+C with columns set as:
#[input,output 1,output 2,alpha,beta,gamma]
#
#xform.dat contains all reactions: A->B with columns set as: 
#[input,output,alpha,beta,gamma]
#
#alpha,beta,gamma represent the Ahrrenius coefficients for the model following UMIST 2012 

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
from sys import exit

#Convert Inputs to Species Number

convert = lambda x: species[str(x,'utf-8')]

#Read Species File 
species={}
with open('data/kida_spec_2016-11-11_1.dat','r') as my_file:
    for line in my_file:
        columns = line.strip().split()
        print(columns)
        species[columns[0]]=int(columns[24])

numspecies=len(species) # Total number of chemical species

print(species,numspecies)

exit()

namesa=("input1","input2","output","alpha","beta","gamma")
namesd=("input1","input2","output1","output2","alpha","beta","gamma")
namesb=("input","output1","output2","alpha","beta","gamma")
namesx=("input","output","alpha","beta","gamma")

assoc    = np.genfromtxt("data/assoc.dat",dtype=None,skip_header=2,names=namesa,converters={0:convert,1:convert,2:convert})
dblassoc = np.genfromtxt("data/dblassoc.dat",dtype=None,skip_header=2,names=namesd,converters={0:convert,1:convert,2:convert,3:convert})
breakup  = np.genfromtxt("data/breakup.dat",dtype=None,skip_header=2,names=namesb,converters={0:convert,1:convert,2:convert})
xform    = np.genfromtxt("data/xform.dat",dtype=None,skip_header=2,names=namesx,converters={0:convert,1:convert}) 

##Arrhenius is the definition which determines the Arrhenius coefficients 
##for each reaction. 

def arrhenius(a,b,c,t):
    a *= 3.1536e7     # Convert units to yr**-1

    def T(t):
        t0 = 0.273790926  # 0.273... years is approximately 100 days
        T0 = 2.1e4        # Initial temperature in Kelvin
        gamma = 1.593     # "adiabatic" index (All data from C&D 2009)
        return T0*(t/t0)**(3-3*gamma)
    
    k=a*(T(t)/300)**b*np.exp(-c/T(t))
    return k

##Defining chemnet 

def chemnet(y,t):
    f=np.zeros(numspecies,float)
##This creates terms with 2 inputs and 1 output
    for a in assoc:
        f[a['input1']] -= arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input1']]*y[a['input2']]
        f[a['input2']] -= arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input1']]*y[a['input2']]
        f[a['output']] += arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input1']]*y[a['input2']]
##This creates terms with 2 inputs and 2 outputs        
    for a in dblassoc:    
        f[a['input1']] -= arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input1']]*y[a['input2']]
        f[a['input2']] -= arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input1']]*y[a['input2']]
        f[a['output1']] += arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input1']]*y[a['input2']]
        f[a['output2']] += arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input1']]*y[a['input2']]        
##This creates terms with 1 input and 2 outputs
    for a in breakup:
        f[a['input']] -= arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input']]
        f[a['output1']] += arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input']]
        f[a['output2']] += arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input']]
##This creates terms with 1 input and 1 output
    for a in xform:
        f[a['input']] -= arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input']]
        f[a['output']] += arrhenius(a['alpha'],a['beta'],a['gamma'],t)*y[a['input']]
    return f
time = np.linspace(1e-10,2,100000)
# Initial Values (Must sum to 1) #########################
yinit = np.zeros(numspecies)          
yinit[0] = 1e5           #Initial Oxygen
yinit[1] = 1e10           #Initial Carbon 
##########################################################
y = odeint(chemnet,yinit,time)

abundance=[]

for i in range(numspecies):
    abundance.append(y[-1,i])

######
# Figure Plotting
######
fig = plt.figure(figsize=(16, 8))
ax1 = fig.add_subplot(121)
ax1.plot(time,y[:,0],label='$O$')
ax1.plot(time,y[:,1],label='$C$')
ax1.plot(time,y[:,2],label='$C_2$')
ax1.plot(time,y[:,3],label='$C_3$')
ax1.plot(time,y[:,4],label='$C_4$')
ax1.plot(time,y[:,5],label='$C_5$')
ax1.plot(time,y[:,6],label='$C_6$')
ax1.plot(time,y[:,7],label='$CO$')
ax1.plot(time,y[:,0]+y[:,1]+2*y[:,2]+3*y[:,3]+4*y[:,4]+5*y[:,5]+6*y[:,6]+2*y[:,7],label='Test')
ax1.set_ylim([0,2.01e10])
ax1.set_xlabel('Time (yr)')
ax1.set_ylabel('Abundance')
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#
ax2=fig.add_subplot(122)
labels = [' ','O','C','C2','C3','C4','C5','C6','CO',' ']
ax2.set_xticklabels(labels)
ax2.plot(abundance,marker="o",linestyle="None")
ax2.set_xlim([-1,numspecies])
ax2.set_yscale('log')
ax2.set_title('Final Abundance')
ax2.set_ylabel('Abundance')
plt.subplots_adjust(left=0.1, right=0.94, top=0.9, bottom=0.1,wspace=0.5)
plt.show()
