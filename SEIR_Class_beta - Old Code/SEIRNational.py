#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import datetime as dt
import Single_dist_ref_SEIR as SDSEIR
import matplotlib.pyplot as plt

#  ---------------------------------- #
#       Find Optimal Parameters       #
# ----------------------------------- #
np.random.seed()
tsim = 400
mov = 0.2
results = SDSEIR.ref_sim_national(mov=mov,qp=0,tsim = tsim,tci=None,movfunct='sawtooth')

Si = results['sim']['S'][0]
Ei = results['sim']['E'][0]
Ii = results['sim']['I'][0]
Ri = results['sim']['R'][0]
ti = 0

beta = results['params'][0]
sigma = results['params'][1]
gamma = results['params'][2]

Ir = pd.DataFrame(results['Ir'],columns=['Ir'])
tr = pd.DataFrame(results['tr'],columns=['tr'])
rd = tr.join(Ir)

params = pd.DataFrame(results['params'].tolist(),columns = ['params'],index=['beta','sigma','gamma','mu'])

print(params)
print('r0 = '+str(beta/gamma))

#aux = pd.read_csv(path+str(i).zfill(2)+'_uni_initdate.csv',index_col=0)   
#initdate = initdate.join(aux)

"""
# ----------------------------- #
#      Quarantine simulation    #
# ----------------------------- #
"""
# No Quarantine
sim = results['sim']['I']
savetocsv=True
# Total Quarantine
qp = -1
aux = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=None, movfunct = 'sawtooth')

sim.append(aux['I'])

sim = pd.DataFrame(sim)

sim = pd.DataFrame(sim)
if savetocsv:
    sim.to_csv('Nacional_sim_TQ.csv')

# Dynamic Quaratine 7D
qp = 7
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=None, movfunct = 'sawtooth')
sim.pop('t')
sim = pd.DataFrame(sim)
if savetocsv:
    sim.to_csv('Nacional_sim_7D.csv')

# Dynamic Quaratine 14D
qp = 14
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=None, movfunct = 'sawtooth')
sim.pop('t')
sim = pd.DataFrame(sim)
if savetocsv:
    sim.to_csv('Nacional_sim_14D.csv')

# Dynamic Quaratine 28D
qp = 28
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=None, movfunct = 'sawtooth')
sim.pop('t')
sim = pd.DataFrame(sim)
if savetocsv:
    sim.to_csv('Nacional_sim_28D.csv')
"""
# ---------------------- #
#      Process Data      #
# ---------------------- #
"""
# Normalize data


"""
# ---------------------- #
#        Plot Data       #
# ---------------------- #
"""
# idx=np.searchsorted(test.t,tr)

# S_dif=S_su1-test.S[:,idx]
# E_dif=E_su1-test.E[:,idx]
# I_dif=I_su1-test.I[:,idx]
# R_dif=R_su1-test.R[:,idx]

# np.amax(S_dif)
# np.amax(E_dif)
# np.amax(I_dif)
# np.amax(R_dif)

# plt.figure()
# plt.plot(test.t[idx],S_dif[1,:],label='Susceptible')
# plt.plot(test.t[idx],E_dif[1,:],label='Exposed')
# plt.plot(test.t[idx],I_dif[1,:],label='Infected simulation')
# plt.plot(test.t[idx],R_dif[1,:],label='Removed')
# plt.xlabel('Days')
# plt.ylabel('Population')
# plt.title('COVID-19 Model')
# plt.legend(loc=0)
# plt.show()

# plt.figure()
# plt.plot(test.t[idx],test.S[1,idx],label='Susceptible')
# plt.plot(test.t[idx],test.E[1,idx],label='Exposed')
# plt.plot(test.t[idx],test.I[1,idx],label='Infected simulation')
# plt.plot(test.t[idx],test.R[1,idx],label='Removed')
# plt.xlabel('Days')
# plt.ylabel('Population')
# plt.title('COVID-19 Model')
# plt.legend(loc=0)
# plt.show()


# plt.figure()
# plt.plot(test.t[idx],S_su1[1,:],label='Susceptible')
# plt.plot(test.t[idx],E_su1[1,:],label='Exposed')
# plt.plot(test.t[idx],I_su1[1,:],label='Infected simulation')
# plt.plot(test.t[idx],R_su1[1,:],label='Removed')
# plt.xlabel('Days')
# plt.ylabel('Population')
# plt.title('COVID-19 Model')
# plt.legend(loc=0)
# plt.show()
"""
# ----------------- #
#    Import Data
# ----------------- #
"""
# NQ, TQ, 7d,14d,28d
data = [[]]

#NQ
sim = pd.read_csv('National_sim_NQ.csv',index_col=0)
dataNQ = sim['I'].to_list()
#TQ
sim = pd.read_csv('Nacional_sim_TQ.csv',index_col=0)
dataTQ = sim['I'].to_list()
#7d
sim = pd.read_csv('Nacional_sim_7D.csv',index_col=0)
data7D = sim['I'].to_list()
#14d
sim = pd.read_csv('Nacional_sim_14D.csv',index_col=0)
data14D = sim['I'].to_list()
#28d
sim = pd.read_csv('Nacional_sim_28D.csv',index_col=0)
data28D = sim['I'].to_list()

data = [dataNQ,dataTQ,data7D,data14D,data28D]

#Real Data
realdata = pd.read_csv('National_RealData.csv')
tci = int(realdata['tr'].to_list()[-1])


params = pd.read_csv('National_params.csv',index_col=0)['params']
beta = params['beta']
sigma = params['sigma']
gamma = params['gamma']
mu = params['mu']

Si = sim['S'][0]
Ei = sim['E'][0]
Ii = sim['I'][0]
Ri = sim['R'][0]
ti = 0

# Simulate

# Total Quarantine
qp = -1
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=tci, movfunct = 'sawtooth')


"""
# ------------------------ #
#   Manage and Save Data   #
# ------------------------ #
"""
savetocsv = True
# Saving real Data
rd.to_csv('National_RealData.csv',index=False)

# Saving parameters
if savetocsv:
    params.to_csv('National_params.csv')

# -------------------------- #
#   Saving Simulation Data   #
# -------------------------- #

# No Quarantine
#sim = results['sim']
#sim.pop('t')
#sim = pd.DataFrame(sim)

qp = 0
tci = 42
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=tci, movfunct = 'sawtooth')
sim.pop('t')
sim = pd.DataFrame(sim)
if savetocsv:
    sim.to_csv('Nacional_sim_NQ.csv')

# Total Quarantine
qp = -1
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=tci, movfunct = 'sawtooth')
sim.pop('t')
sim = pd.DataFrame(sim)
sim.to_csv('Nacional_sim_TQ.csv')

# Dynamic Quaratine 7D
qp = 7
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=tci, movfunct = 'sawtooth')
sim.pop('t')
sim = pd.DataFrame(sim)
if savetocsv:
    sim.to_csv('Nacional_sim_7D.csv')

# Dynamic Quaratine 14D
qp = 14
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=tci, movfunct = 'sawtooth')
sim.pop('t')
sim = pd.DataFrame(sim)
if savetocsv:
    sim.to_csv('Nacional_sim_14D.csv')

# Dynamic Quaratine 28D
qp = 28
sim = SDSEIR.intger(Si,Ei,Ii,Ri,ti,tsim,0.01,beta,sigma,gamma,mov,qp=qp,tci=tci, movfunct = 'sawtooth')
sim.pop('t')
sim = pd.DataFrame(sim)
if savetocsv:
    sim.to_csv('Nacional_sim_28D.csv')





""" 
-----------------------------------------------------------------------------------------------------------------------
Pipeline Nacional oficial
"""

realdata = pd.read_csv('TotalesNacionales.csv')
Ireal = realdata.loc[4].to_list()[1:]
