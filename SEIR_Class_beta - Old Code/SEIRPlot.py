#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Single_dist_ref_SEIR as SDSEIR


"""
 
 Plot pipeline for SEIR simulations


"""
# Select state and comuna: 
state = '05'
comuna = '05101'

# ---------------------------- # 
#   Import and process Data    #
# ---------------------------- #

# Paths
path = '../Data/DatosConsolidados/'
filename_params = '05_uni_parameters.csv'
filename_data = 'ActivosAcumuladosAle.xlsx'
filename_cutcomuna = 'cut_2018_v03.xls'

# Import parameters
parameters = pd.read_csv(path+filename_params,index_col=0)
beta = parameters[comuna]['beta']
gamma = parameters[comuna]['gamma']
sigma = parameters[comuna]['sigma']
mu = parameters[comuna]['mu']
err = parameters[comuna]['err']

# Import real data
data = pd.read_excel(path+filename_data,index_col=0).transpose()

# --------------------------- #
#   Date number calculation   #
# ----------------------------#
a = data.iterrows()
aux = [i[0] for i in a]
tr = [(aux[i]-aux[0]).days for i in range(len(aux))]

# ----------------------- #
#   Cut to Comuna Name    #
# ----------------------- #
# to do mejorar la busqueda para las comunas en que los nombres no calcen 
datacut = pd.read_excel(path+filename_cutcomuna)
try:
    ncomuna = datacut.loc[datacut['Código Comuna 2017']==int(comuna),'Nombre Comuna'].tolist()[0]
    Ir = data[ncomuna].tolist()
except:
    print('Comuna '+ncomuna+' not found in data')
    raise Exception('Comuna name not found')


# ------------------------------- #
#   From accumulated to actives   #
# ------------------------------- #

# New Cases
Cn=np.zeros(len(data[ncomuna]))
Cn[0]=data[ncomuna][0]
for i in range(1,len(Cn)):
    Cn[i]=data[ncomuna][i]-data[ncomuna][i-1]

# Total active cases each day
Ir=np.zeros(data.shape[0])
Ir[np.where(np.array(tr)<=14)[0]]=data[ncomuna][np.where(np.array(tr)<=14)[0]]
for i in np.where(np.array(tr)>14)[0]:
    ind=np.where(np.array(tr)-tr[i]+14<0)[0]
    Ir[i]=data[ncomuna][i]-sum(Cn[0:ind[-1]])    



qp = 0
mov = 0.2 
tsim = 100
tci = None
movfunct='sawtooth'


# refine again madafacas
refresults = SDSEIR.ref_sim_cons(state,comuna,mov,qp,tsim,tci,movfunct)
beta = refresults['params'][0]
gamma = refresults['params'][1]
sigma = refresults['params'][2]
mu = refresults['params'][3]
err = refresults['err']

# ------------ #
#   Simulate   #
# ------------ #

results = SDSEIR.simulate(state,comuna,beta,sigma,gamma,mu,qp,mov,tsim,tci,movfunct)
results['init_date'] = '22-feb-2020'
S = results['S']
E = results['E']
I = results['I']
R = results['R']
t = results['t']
# results keys = ['S', 'E', 'I', 'R', 't', 'init_date'])


# ---------- #
#    Plot    #
# ---------- #

plt.figure()
#plt.plot(t,S,label='Susceptible')
# plt.plot(t,E,label='Exposed')
plt.plot(t,I,label='Infected simulation')
plt.plot(tr,Ir,label='Infected real',color ='red',marker='o')
# plt.plot(t,R,label='Removed')
# plt.plot(t,R,label='Removed')
plt.xlabel('Days')
plt.ylabel('Population')
plt.title('COVID-19 Model '+ncomuna)
plt.legend(loc=0)
plt.show()


