# -*- coding: utf-8 -*-
"""
input.py - Input Reaction Data into CarBoN

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

"""

import pandas as pd

"""
This code inputs a KIDA reactions and species file using Pandas. Column names are provided by KIDA.

Outputs:  

reactions - dataframe for all reactions in the network
speciesidx - dictionary of species name and index
speciesmass - dictionary of species name and mass*
numspecies - Total number of chemical species + a "catchall" zero.


*Mass is calculated per atom, i.e. O2 has a mass of 2 and Si2O2 has a mass of 4. 

"""
def kida_input(reac_file,spec_file):
    speciesidx={}
    speciesmass={}
    colnames=["Input1","Input2","Input3","Output1","Output2","Output3","","","alpha","beta","gamma","F","g",\
                  "Type_of_uncertainty","itype","T_low","T_high","Formula","Number","Numbers of","Recommendation"]
    colwidths=[11,11,12,11,11,11,11,12,11,11,11,9,9,5,3,7,7,3,5,2,3]
    collist=[0,1,3,4,5,8,9,10,11,12,13,14,15,16,17,18,19,20]
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


#Convert Inputs to Species Number
#Some inputs may be missing (typically Output2 or Output3) Don't change those values. Let Pandas do that.

    convert = lambda x: int(speciesidx[x]) if x!='' else None


#Read the reactions file and replace species names with index numbers from species file.
#First two lines read the file with properly formatted columns. 
#See: http://stackoverflow.com/questions/40663510/pandas-read-fwf-ignoring-values
#
#Force index number columns to be integers. reactions.fillna replaces NaN with float.

    reactions=pd.read_fwf(reac_file,comment='#',names=colnames,widths=colwidths,usecols=collist,converters={0:convert,1:convert,2:convert,3:convert,4:convert})
    reactions.fillna(0,inplace=True)
    reactions[["Input1","Input2","Output1","Output2","Output3"]] = reactions[["Input1","Input2","Output1","Output2","Output3"]].astype(int)
    
    return reactions,speciesidx,speciesmass,numspecies

"""
This function prints the entire Pandas dataframe (to test quality of inputs)

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
