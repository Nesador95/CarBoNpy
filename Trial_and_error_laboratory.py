# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 11:32:21 2019

@author: Luis Felipe Llano (Python 3.6)

Desciption:
"""
import numpy as np
import matplotlib.pyplot as plt
import datainput






datainput.KIDA_spec(r'data/kida_spec_C_O_only.dat')
df = datainput.KIDA_spec(r'data/kida_spec_C_O_only.dat')
df1 = datainput.KIDA_reac(r'data/kida_reac_C_O_only.dat')

final_table = datainput.KIDA_reac_com_readable(df1, df)
#df['pivot_columns'] = [0 for i in list(df.index)]


#print(df)
#print(list(df.index))
#print(df1)
#print(len(df.index))
#print(np.zeros([len(df.index)]))
#print(np.zeros([len(df.index),len(df.columns)]))
#print(len(df.index))

#series = df[['species','species_#','pivot_columns']]
#series = series.pivot(columns= 'pivot_columns',
#             index='species',
#             values='species_#')
#print(series)
#series1 = series[0]
#print(type(series1))
#print(series1)
#series1 = series.to_dict()
#print(series1)
#dictionary = series1[0]
#print(dictionary)


#df1['Input1_id'] = df1['Input1'].map(dictionary)
#print(df1)
#print(df[['species','species_id']])

#print(final_table)
#simple_final_table = final_table[['Input1', 'Input1_id','Input2', 'Input2_id',
#                   'Output1', 'Output1_id','Output2', 'Output2_id',
#                   'Output3', 'Output3_id']]
#print(simple_final_table)
#def counter(data):
#    count = 0
#    for i in list(data):
#        count += 1
#    return count

#print(counter(kida_spec.index))
#print(list(kida_spec['species_#']))
#print(len(kida_spec.index))
#species = 2
#print(kida_spec['species_#'].values[species])

print(len(list(df.index)))
count=len(list(df.index))
print(count)
print(range(count))
for num in range(count):
        
    in1=final_table.loc[num]['Input1_id']
    in2=final_table.loc[num]['Input2_id']
    #in3=int(kida_reac.loc[num,'Input3_id'])
    out1=final_table.loc[num]['Output1_id']
    out2=final_table.loc[num]['Output2_id']
    out3=final_table.loc[num]['Output3_id']
    alpha=final_table.loc[num]['alpha']
    beta=final_table.loc[num]['beta']
    gamma=final_table.loc[num]['gamma']
    fo=final_table.loc[num]['Formula']
    print(in1)

print(df['species_#'].to_dict())               