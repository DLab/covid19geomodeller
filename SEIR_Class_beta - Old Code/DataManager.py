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


path = '../Data/unirefine/'

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

simdata.to_csv(path+'sim_02.csv')
parameters.to_csv(path+'parameters_02.csv')
initdate.to_csv(path+'initdate_02.csv')

"""
# Correjir fecha:

cutlist = []
cutlistpath = "../Data/cutlist.csv"
cutlist = pd.read_csv(cutlistpath, header = None,dtype=str)

initdate = pd.DataFrame()

# Refine parameters per CUT    

for index, row in cutlist.iterrows():    
    state = str(row[0])[0:2]
    comuna = str(row[0])




endpoint="http://192.168.2.220:8080/covid19/findComunaByIdState?idState="+state+"&&comuna="+comuna
r = requests.get(endpoint) #, params = {"w":"774508"})
mydict = r.json()
info=pd.DataFrame(mydict)

endpoint="http://192.168.2.220:8080/covid19/getDatosMinSalSummary?state="+state+"&comuna="+comuna
r = requests.get(endpoint) #, params = {"w":"774508"})
mydict = r.json()
data=pd.DataFrame(mydict)
data=data[data.data != 0]
data=data.reset_index()
"""