#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 12:55:18 2020

@author: sropert

This script uses the OD matrix to locate the different urban centers that are isolated between each other 
"""
import pandas as pd
import numpy as np

# ---------------------- #
#    Import OD Matrix    #
# ---------------------- #
path = "../Data/UrbanCenters/"
odfile = 'Movilidad_diaria_comuna.csv'
OD = pd.read_csv(path+odfile, index_col=0)

comunas = list(OD.keys())
urbancenters = dict()
j = 0
for i in comunas:
    aux = OD[i] == OD[i]
    # Districts connected to district i
    if True in aux.tolist():
        print('True in aux')
        urbancenters[j] = [k for k in comunas if aux[k]]
        #urbancenters[j].append(i)
        for l in urbancenters[j]:
            aux2 = OD[l] == OD[l]
            for k in comunas:
                if aux2[k] and (k not in urbancenters[j]):
                    urbancenters[j].append(k)

        # Remove from comunas
        try:
            for k in urbancenters[j]:
                comunas.remove(k)
        except:
            pass
        j+=1
    else:
        comunas.remove(i)
    

# Create new OD matrix for each urban center
urbancentersDF = dict()
for i in urbancenters: 
    urbancentersDF[i] = OD[urbancenters[i]].loc[urbancenters[i]]
    filename = path+'UrbanCenter_'+str(i)+'.csv'
    urbancentersDF[i].to_csv(filename)
