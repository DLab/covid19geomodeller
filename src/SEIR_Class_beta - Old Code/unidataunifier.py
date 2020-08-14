#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from class_SEIR import SEIR
from  SEIRrefiner import SEIRrefiner 
import pandas as pd
import numpy as np
from timeit import default_timer as timer
import Single_dist_ref_SEIR as SDSEIR
import csv
import datetime


#path = '../Data/unirefine/'
path = '../Data/DatosConsolidados/'
# import data list
# init dataframe
simdata = pd.read_csv(path+'01_uni_sim.csv',index_col=0)
parameters = pd.read_csv(path+'01_uni_parameters.csv',index_col=0)
initdate = pd.read_csv(path+'01_uni_initdate.csv',index_col=0)

for i in range(2,17):
    aux = pd.read_csv(path+str(i).zfill(2)+'_uni_sim.csv',index_col=0)
    simdata = simdata.join(aux)
    aux = pd.read_csv(path+str(i).zfill(2)+'_uni_parameters.csv',index_col=0)
    parameters = parameters.join(aux)
    aux = pd.read_csv(path+str(i).zfill(2)+'_uni_initdate.csv',index_col=0)   
    initdate = initdate.join(aux)

simdata.to_csv(path+'sim.csv')
parameters.to_csv(path+'parameters.csv')
initdate.to_csv(path+'initdate.csv')

#pd.read_csv("/media/samuel/Data/samue/Dropbox/DLab/Data/Unicomunal/parameters.csv", index_col=0)   
