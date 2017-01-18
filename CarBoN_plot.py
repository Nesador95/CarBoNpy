#!/usr/bin/python
#title			:CarBon_plot.py
#description		:Plots Saved CarBoN Data
#author			:Ethan Deneault
#date			:20161220
#version		:0.1
#usage			:python CarBoN_plot.py
#notes			:
#python_version		:3.5.1
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from sys import exit

outfile = 'output/D_2017/CD_Model/output_test.dat'

infile = np.load(outfile+".npz")

print(infile.files)

y=infile['y']/1.1
speciesidx = infile['speciesidx'][()]  #extracts dictionary from array?
speciesmass =infile['speciesmass'][()]
abundance=infile['abundance']
time=365.25*(infile['time']-0.1743)

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


######
# Figure Plotting
######
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111)
jet = cm = plt.get_cmap('nipy_spectral') 
cNorm  = colors.Normalize(vmin=0, vmax=8)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

def colr(n):
    stop = False
    if n in [0,1,2,3,4,5,6]:
        return n
    while stop != True:
        n -= 6
        if n in [1,2,3,4,5,6]:
            stop = True
            return n
        
for name,index in speciesidx.items():
#    linecolor=colr(index)
#    colorVal = scalarMap.to_rgba(linecolor)
    if index in [1,2]:
        ax1.plot(time,y[:,index],label=name,linewidth=2.0)
    if index in [3,4]:
        ax1.plot(time,y[:,index],'--',label=name,linewidth=2.0)
    if index in [5,6]:
        ax1.plot(time,y[:,index],':',label=name,linewidth=2.0)
#    if index in [17,18,19,20]:
#        ax1.plot(time,y[:,index],ls='dashdot',label=name,linewidth=2.0) 

plt.subplots_adjust(left=0.1, right=0.84, top=0.9, bottom=0.1,wspace=0.5) 
ax1.set_ylim([1e-12,3e0])
ax1.plot(time,total,label='Test')
#ax1.set_xlim([0,1000])
ax1.set_xlim([0,1760])
ax1.set_yscale('log')
ax1.set_xlabel('t-t0 (days)')
ax1.set_ylabel('Abundance')
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

####################################################

#ax1.plot(time,total,label='Test')
ax1.set_ylim([5e-11,2.01e10])
#ax1.set_xlim([73,1826])
#ax1.set_xlim([550,600])
ax1.set_yscale('log')
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Abundance')
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
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
plt.show()
