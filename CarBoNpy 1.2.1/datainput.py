# -*- coding: utf-8 -*-
"""
input.py - Input and print subroutines for CarBoN

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

"""

import pandas as pd
import configparser as cp

"""
This function reads in a KIDA reactions and species file using Pandas.
Column names and widths are provided by KIDA.

Outputs:  

reactions - dataframe for all reactions in the network
speciesidx - dictionary of species name and index
speciesmass - dictionary of species name and mass*
numspecies - Total number of chemical species + a "catchall" zero.


*Mass is calculated per atom, i.e. O2 has a mass of 2
and Si2O2 has a mass of 4. 

Big Thanks to: 
http://stackoverflow.com/questions/40663510/pandas-read-fwf-ignoring-values

"""
##############################################################################
# The following loads the imput from a KIDA file
##############################################################################

def KIDA_reac(reac_file):
    
    colwidths=[11,23,11,10,35,10,11,11,9,9,6,3,8,6,4,4,3,2]

    reac_df = pd.read_fwf(reac_file, comment='#', widths=colwidths)
    reac_df.rename(columns={"Type": "Type of uncertainty",
                       "Re": "itype",
                       "Tlo": "T low",
                       "Thi": "T high",
                       "Fo": "Formula",
                       "N": "Number",
                       "V": "Number of",
                       "R": "Recomendation"}, inplace=True)
    return reac_df


def KIDA_spec(spec_file):
    
    colwidths=[11,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5]
    
    spec_df = pd.read_fwf(spec_file, comment='#', widths=colwidths, 
                          header=None)
    spec_df.columns = list(range(0,len(colwidths)))
    
    spec_df['#_of_atoms'] = spec_df[2]+spec_df[3]+spec_df[4]+\
            spec_df[5]+spec_df[6]+spec_df[7]+spec_df[8]+spec_df[9]+\
            spec_df[10]+spec_df[11]+spec_df[12]+spec_df[13]+spec_df[14]\
            +spec_df[15]+spec_df[16]+spec_df[17]+spec_df[18]+spec_df[19]\
            +spec_df[20]+spec_df[21]+spec_df[22]+spec_df[23]
    spec_df = spec_df.drop([2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                            20,21,22,23], axis=1)
    spec_df.rename(columns={0 : "species",
                       1 : "charge",
                       24: "species_#"}, inplace=True)
    
    spec_df.at[18, 'species_#'] = 19
    spec_df.at[19, 'species_#'] = 20
    spec_df.at[20, 'species_#'] = 21                     
    spec_df.at[21, 'species_#'] = 99
               
    spec_df['species_#'] = spec_df['species_#'] - 1
    
                                
    return spec_df

def KIDA_reac_com_readable(KIDA_reac, KIDA_spec):
    KIDA_spec['pivot_column'] = [0 for i in list(KIDA_spec.index)]
    
    series = KIDA_spec[['species','species_#','pivot_column']]
    series = series.pivot(columns= 'pivot_column',
             index='species',
             values='species_#')
    series1 = series[0]
    series1 = series.to_dict()
    
    dictionary = series1[0]
    
    KIDA_reac['Input1_id'] = KIDA_reac['Input1'].map(dictionary)
    KIDA_reac['Input2_id'] = KIDA_reac['Input2'].map(dictionary)
    KIDA_reac['Output1_id'] = KIDA_reac['Output1'].map(dictionary)
    KIDA_reac['Output2_id'] = KIDA_reac['Output2'].map(dictionary)
    KIDA_reac['Output3_id'] = KIDA_reac['Output3'].map(dictionary)
    
    KIDA_reac['Input1_id'] = KIDA_reac['Input1_id'].fillna(0)
    KIDA_reac['Input2_id'] = KIDA_reac['Input2_id'].fillna(0)
    KIDA_reac['Output1_id'] = KIDA_reac['Output1_id'].fillna(0)
    KIDA_reac['Output2_id'] = KIDA_reac['Output2_id'].fillna(0)
    KIDA_reac['Output3_id'] = KIDA_reac['Output3_id'].fillna(0)
    
    return KIDA_reac
##############################################################################
# The following loads the imput from a BASIC file
##############################################################################

def basic_input(reac_file,spec_file):
    """
    This function reads in "basic format" reactions/species files using Pandas. 
    Basic Reaction Format is a comma delimited file. 

    Input1 Input2 Input3 Output1 Output2 Output3 Alpha Beta Gamma Formula

    For all Arrhenius reactions, Formula should be 3.

    Basic Species Format:

    Name Mass Index

    This probably needs more documentation. 

    """
    speciesidx={}
    speciesmass={}
    colnames=["Input1", "Input2", "Input3", "Output1", "Output2", "Output3",
              "alpha", "beta", "gamma", "Formula"]
    with open(spec_file,'r') as my_file:
        for line in my_file:
            columns = line.strip().split()
            speciesidx[columns[0]]=int(columns[2])
            speciesmass[columns[0]]=int(columns[1])
                
    numspecies=len(speciesidx)+1

    convert = lambda x: speciesidx[x] if x!=' ' else None

    reactions=pd.read_csv(reac_file, comment = '#', names = colnames,
                          converters = {0:convert, 1:convert, 2:convert,
                                        3:convert, 4:convert, 5:convert})
    reactions.fillna(0,inplace = True)
    reactions[colnames[0:6]] = reactions[colnames[0:6]].astype(int)
    
    return reactions,speciesidx,speciesmass,numspecies

##############################################################################
# The following prints the entire list straight from pandas
##############################################################################

def print_full(x):
    """
    This function prints the entire Pandas dataframe to test. 
    Make sure that your terminal window is wide enough!

    Code shamlessly absorbed from stackoverflow
    http://stackoverflow.com/questions/19124601/is-there-a-way-to-pretty-print-the-entire-pandas-series-dataframe

    """
    pd.set_option('display.max_rows', len(x))
    pd.set_option('display.max_columns', 1000)
    pd.set_option('display.expand_frame_repr', False)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.expand_frame_repr')

##############################################################################
# The following reads the settings file
##############################################################################

def settings():
    
    """ This function reads the settings specified in settings.ini """

    config = cp.ConfigParser()
    config.read('settings.ini')
    
    files = config['files']
    model = config['model']
    plot = config['plot']
    
    file_format = files['format']
    species_file = r'data/' + files['species file']
    reactions_file = r'data/' + files['reactions file']
    output_file = 'output/' + files['output file']
    model_type = model['model type']
    density = model.getfloat('density')
    temperature = model.getfloat('temperature')
    start_time = model.getfloat('start time')
    end_time = model.getfloat('end time')
    outfile = plot['outfile for plotting']

    return file_format, species_file, reactions_file, output_file, model_type, density, temperature, start_time, end_time, outfile

##############################################################################
# The following reads the abundances file
##############################################################################

def abundances():

    """ This function reads the abundances.ini file """
    
    initial = {}
    with open('abundances.ini','r') as my_file:
        header = my_file.readline()
        for line in my_file:
            columns = line.strip().split()
            initial[columns[0]] = float(columns[1])
            
    return initial
