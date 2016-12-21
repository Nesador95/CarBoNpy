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
from sys import exit

outfile = 'output/kida_output_C_O_only.dat'

infile = np.load(outfile+".npz")

print(infile.files)

y=infile['y']
speciesidx = infile['speciesidx'][()]  #extracts dictionary from array?
abundance=infile['abundance']
time=365.25*infile['time']

total = y[:,1]+y[:,2]+2*y[:,3]+2*y[:,4]+3*y[:,5]+4*y[:,6]+5*y[:,7]+6*y[:,8]+7*y[:,9]+8*y[:,10]+9*y[:,11]+10*y[:,12]+2*y[:,13]

print(speciesidx)

######
# Figure Plotting
######
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(111)
for name,index in speciesidx.items():
    if index in [1,2,3,4,5]:
        ax1.plot(time,y[:,index],label=name,linewidth=2.0)
#    if index in [5,6,7,8]:
#        ax1.plot(time,y[:,index],'--',label=name,linewidth=2.0)
#    if index in [9,10,11,12]:
#       ax1.plot(time,y[:,index],':',label=name,linewidth=2.0)
#    if index in [13,14,15,16]:
#        ax1.plot(time,y[:,index],ls='dashdot',label=name,linewidth=2.0) 

plt.subplots_adjust(left=0.1, right=0.84, top=0.9, bottom=0.1,wspace=0.5) 
ax1.set_ylim([1e-12,2.01e10])
ax1.plot(time,total,label='Test')
ax1.set_xlim([73,1826])
#ax1.set_xlim([550,600])
ax1.set_yscale('log')
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Abundance')
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

####################################################

#ax1.plot(time,total,label='Test')
ax1.set_ylim([1e-12,2.01e10])
ax1.set_xlim([73,1826])
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
