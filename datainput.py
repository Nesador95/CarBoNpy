# -*- coding: utf-8 -*-
"""
input.py - Input and print subroutines for CarBoN

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

"""

import pandas as pd
import configparser as cp

"""
This function reads in a KIDA reactions and species file using Pandas. Column names and widths are provided by KIDA.

Outputs:  

reactions - dataframe for all reactions in the network
speciesidx - dictionary of species name and index
speciesmass - dictionary of species name and mass*
numspecies - Total number of chemical species + a "catchall" zero.


*Mass is calculated per atom, i.e. O2 has a mass of 2 and Si2O2 has a mass of 4. 

Big Thanks to: http://stackoverflow.com/questions/40663510/pandas-read-fwf-ignoring-values
"""

def kida_input(reac_file,spec_file):
    speciesidx={}
    speciesmass={}
    colnames=["Input1","Input2","Input3","Output1","Output2","Output3","","","alpha","beta","gamma","F","g",\
                  "Type_of_uncertainty","itype","T_low","T_high","Formula","Number","Numbers of","Recommendation"]
    colwidths=[11,11,12,11,11,11,11,12,11,11,11,9,9,5,3,7,7,3,5,2,3]
    collist=[0,1,2,3,4,5,8,9,10,11,12,13,14,15,16,17,18,19,20]
    with open(spec_file,'r') as my_file:
        for line in my_file:
            columns = line.strip().split()
            speciesidx[columns[0]]=int(columns[24])
            sum=0
            for i in range(2,24):
                num = int(columns[i])
                sum += num
            speciesmass[columns[0]]=sum 
     
    numspecies=len(speciesidx)+1

    convert = lambda x: speciesidx[x] if x!='' else None

    reactions=pd.read_fwf(reac_file,comment='#',names=colnames,widths=colwidths,usecols=collist,\
                              converters={0:convert,1:convert,2:convert,3:convert,4:convert,5:convert})
    reactions.fillna(0,inplace=True)
    reactions[colnames[0:6]] = reactions[colnames[0:6]].astype(int)
    
    return reactions,speciesidx,speciesmass,numspecies

"""
This function reads in "basic format" reactions/species files using Pandas. 
Basic Reaction Format is a comma delimited file. 

Input1 Input2 Input3 Output1 Output2 Output3 Alpha Beta Gamma Formula

For all Arrhenius reactions, Formula should be 3.

Basic Species Format:

Name Mass Index

This probably needs more documentation. 
"""

def basic_input(reac_file,spec_file):
    speciesidx={}
    speciesmass={}
    colnames=["Input1","Input2","Input3","Output1","Output2","Output3","alpha","beta","gamma","Formula"]
    with open(spec_file,'r') as my_file:
        for line in my_file:
            columns = line.strip().split()
            speciesidx[columns[0]]=int(columns[2])
            speciesmass[columns[0]]=int(columns[1])
                
    numspecies=len(speciesidx)+1

    convert = lambda x: speciesidx[x] if x!=' ' else None

    reactions=pd.read_csv(reac_file,comment='#',names=colnames,converters={0:convert,1:convert,2:convert,3:convert,4:convert,5:convert})
    reactions.fillna(0,inplace=True)
    reactions[colnames[0:6]] = reactions[colnames[0:6]].astype(int)
    
    return reactions,speciesidx,speciesmass,numspecies

"""
This function prints the entire Pandas dataframe to test. Make sure that your terminal window is wide enough!

Code shamlessly absorbed from stackoverflow
http://stackoverflow.com/questions/19124601/is-there-a-way-to-pretty-print-the-entire-pandas-series-dataframe
"""

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    pd.set_option('display.max_columns', 1000)
    pd.set_option('display.expand_frame_repr', False)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.expand_frame_repr')

"""
This function reads the settings.ini file. 
"""

def settings():
    config = cp.ConfigParser()
    config.read('settings.ini')
    
    files=config['files']
    model=config['model']
    
    file_format = files['format']
    species_file = 'data/' + files['species file']
    reactions_file = 'data/' + files['reactions file']
    output_file = 'output/' + files['output file']
    model_type = model['model type']
    density = model.getfloat('density')

    return file_format,species_file,reactions_file,output_file,model_type,density

"""
This function reads the abundances.ini file
"""

def abundances():
    initial = {}
    with open('abundances.ini','r') as my_file:
        header = my_file.readline()
        for line in my_file:
            columns = line.strip().split()
            initial[columns[0]]=float(columns[1])
            
    return initial
