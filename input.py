# -*- coding: utf-8 -*-
"""
input.py - Input Reaction Data into CarBoN

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

"""

import pandas as pd

def kida_input(reac_file,spec_file):
    speciesidx={}
    speciesmass={}
    colnames=["Input1","Input2","Output1","Output2","Output3","alpha","beta","gamma","F","g",\
                  "Type_of_uncertainty","itype","T_low","T_high","Formula","Number","Numbers of","Recommendation"]
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

    convert = lambda x: int(speciesidx[x]) if x!='' else None

    reactions=pd.read_fwf(reac_file,comment='#',names=colnames,converters={0:convert,1:convert,2:convert,3:convert,4:convert})
    reactions.fillna(0,inplace=True)
    reactions[["Input1","Input2","Output1","Output2","Output3"]] = reactions[["Input1","Input2","Output1","Output2","Output3"]].astype(int)
    return reactions,speciesidx,speciesmass,numspecies
