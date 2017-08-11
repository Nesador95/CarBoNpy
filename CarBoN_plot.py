#!/usr/bin/python
#title			:CarBon_plot.py
#description		:Plots Saved CarBoN Data
#author			:Ethan Deneault1
#date			:20161220
#version		:0.1
#usage			:python CarBoN_plot.py
#notes			:
#python_version		:3.5.1
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#import matplotlib._color_data 
#import matplotlib.cm as cmx
#from sys import exit

outfile = 'basic_C-0.01_O-1_Si-0.1_Compton_D1e10_T1.8e4.dat'
infile = np.load(outfile+".npz")

totmass = 1.11 ### Total of C+O+Si

plottitle = 'C/O = 1/100      C/Si = 1/10'

print(infile.files)

y=infile['y']/totmass
speciesidx = infile['speciesidx'][()]  #extracts dictionary from array?
speciesmass =infile['speciesmass'][()]
abundance=infile['abundance']
time=365.25*(infile['time'])#-(99/365.25))

#print(speciesidx)
#print(y[-1,11])

######
# Tests for mass conservation. Total must = 1. 
######

total=0
for name,index in speciesidx.items():
    if index!=99:
        total += speciesmass[name]*y[:,index]

#exit()

colors=['xkcd:red','xkcd:blue','xkcd:green','xkcd:cyan','xkcd:magenta',\
'xkcd:goldenrod','xkcd:orange','xkcd:indigo']

######
# Figure Plotting
######
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111)

ax1.xkcd()
ax1.plot(time,y[:,4],label='$\mathrm{CO}$',linewidth=2.0,color=colors[0])
ax1.plot(time,y[:,5],label='$\mathrm{C_2}$',linewidth=2.0,color=colors[1])
ax1.plot(time,y[:,6],label='$\mathrm{C_3}$',linewidth=2.0,color=colors[2])
ax1.plot(time,y[:,14],label='$\mathrm{O_2}$',linewidth=2.0,color=colors[3])
ax1.plot(time,y[:,12],'--',label='$\mathrm{SiC}$',linewidth=2.0,color=colors[4])
ax1.plot(time,y[:,13],'--',label='$\mathrm{SiO}$',linewidth=2.0,color=colors[5])
ax1.plot(time,y[:,15],'-.',label='$\mathrm{(SiC)_2}$',linewidth=2.0,color=colors[6])
ax1.plot(time,y[:,16],'-.',label='$\mathrm{(SiO)_2}$',linewidth=2.0,color=colors[7])

#ax1.plot(time,total,label='Test',color='black')
        
#for name,index in speciesidx.items():
#    if index in [4,5,6,7]:
#        ax1.plot(time,y[:,index],label=name,linewidth=2.0,color=colors[index-4])
#    if index in [12,13,14,15,16]:
#        ax1.plot(time,y[:,index],'--',label=name,linewidth=2.0,color=colors[index-12])
#    if index in [12]:
#        ax1.plot(time,y[:,index],':',label=name,linewidth=2.0,color=colors[index-12])
#    if index in [17,18,19,20]:
#        ax1.plot(time,y[:,index],ls='dashdot',label=name,linewidth=2.0) 

plt.subplots_adjust(left=0.1, right=0.84, top=0.9, bottom=0.1,wspace=0.5) 
ax1.set_ylim([1e-20,2e-2])
#ax1.set_xlim([0,200])
ax1.set_xlim([100,1760])
ax1.set_yscale('log')
ax1.set_xlabel('t-t0 (days)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax1.set_ylabel('Abundance')
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax1.set_title(plottitle)
plt.show()

####################################################

#ax1.plot(time,total,label='Test')
#ax1.set_ylim([5e-11,2.01e10])
#ax1.set_xlim([73,1826])
#ax1.set_xlim([550,600])
#ax1.set_yscale('log')
#ax1.set_xlabel('Time (days)')
#ax1.set_ylabel('Abundance')
#ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#
#ax2=fig.add_subplot(122)
#labels = [' ','O','C','C2','C3','C4','C5','C6','CO',' ']
#ax2.set_xticklabels(labels)
#ax2.plot(abundance,marker="o",linestyle="None")
#ax2.set_xlim([-1,numspecies])
#ax2.set_yscale('log')
#ax2.set_title('Final Abundance')
#ax2.set_ylabel('Abundance')
#plt.subplots_adjust(left=0.1, right=0.94, top=0.9, bottom=0.1,wspace=0.5)
#plt.show()
