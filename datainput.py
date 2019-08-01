# -*- coding: utf-8 -*-
"""
input.py - Input and print subroutines for CarBoN

Copyright (c) 2012-2016 Ethan Deneault

This file is part of CarBoN

"""

import pandas as pd
import configparser as cp
from sys import exit


def KIDA_input(reac_file,spec_file):
    '''
    Read in a KIDA (https://kida.obs.u-bordeaux1.fr) formatted reactions 
    and species files.

    Each species in the system is indexed by number, and the output dataframe 
    replaces the species name with the appropriate number.

    Inputs: KIDA-supplied reactions file
            KIDA-supplied species file
    
    Outputs: Dataframe for reactions
             Dataframe for species -> dictionary for assigning species to 
                                      a row in the solution matrix. 
 
    '''
    reac_col_names=['Input1','Input2','Output1','Output2','Output3',
                    'alpha','beta','gamma','F','g','Type','Re','Tlo', 
                    'Thi','Fo','N','V','R']
    reac_col_widths=[11,23,11,11,34,11,11,11,9,9,5,3,7,7,3,5,2,3]
    reac_dtypes={'Input1':str,'Input2':str,'Output1':str,'Output2':str,
                 'Output3':str}

    reac_df = pd.read_fwf(reac_file, comment='#', header=None,
                          names=reac_col_names,widths=reac_col_widths,
                          converters=reac_dtypes)

    spec_df = pd.read_fwf(spec_file, comment='#', header=None)

    col_list=list(spec_df)

    spec_df['atom_num'] = spec_df[col_list[2:24]].sum(axis=1)
    spec_df = spec_df.drop(col_list[2:24], axis=1)
    spec_df.rename(columns={0 : "species",
                       1 : "charge",
                       24: "species_num"}, inplace=True)
    #spec_df.loc[spec_df.species == 'e-', 'species_num'] = 0

    spec_dict = pd.Series(spec_df.species_num.values,
                           index=spec_df.species).to_dict()
    add_photons={'Photon':0,'Pho':0}
    spec_dict.update(add_photons)

    print(spec_dict)

    for column in reac_df.columns[0:5]: 
        reac_df[column] = reac_df[column].replace(spec_dict)
        reac_df[column] = reac_df[column].astype('Int64')

    print(reac_df.dtypes)

    return reac_df, spec_df, spec_dict


def settings():
    
    """ 
    This function reads the settings specified in settings.ini 
    """

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

def abundances(dictionary):

    """ This function reads the initial abundances of reactants supplied 
    by the abundances.ini file 
    """
    
    abund_df = pd.read_fwf('abundances.ini', comment='#')
    abund_df['Species'] = abund_df['Species'].map(dictionary)
            
    return abund_df


def basic_input(reac_file,spec_file):
    """
    This function reads in "basic format" reactions/species files using Pandas. 
    Basic Reaction Format is a comma delimited file. 

    Input1 Input2 Input3 Output1 Output2 Output3 Alpha Beta Gamma Formula

    For all Arrhenius reactions, Formula should be 3.

    Basic Species Format:

    Name Mass Index

    This probably needs more documentation, and can likely be depreciated?

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


