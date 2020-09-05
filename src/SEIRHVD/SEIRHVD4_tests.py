#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

import sys
from pathlib import Path
sys.path.insert(1, '../src/SEIRHVD/')
sys.path.insert(1, '../src/utils/')
sys.path.insert(1, 'src/SEIRHVD/')
sys.path.insert(1, 'src/utils/')

from class_SEIRHUVD4 import SEIRHVD 
from Quarantine import Quarantine
from importdata import ImportData

"""
SEIRHVD 4 Tests & QA

"""

tsim = 1000
beta = 0.2
mu = 1.5

k=0
Htot=30
Vtot=20
H0=0
V0=0
B0=0
D0=0
R0=0
I0=100
I_d0=10
I_ac0=100
SeroPrevFactor=1
expinfection=0
population=1000000

# Quarantines 
max_mob = 0.8
rem_mob = 0.5

alpha = Quarantine(rem_mob).alpha


simulation = SEIRHVD(tsim,beta,mu,alpha,k=k,Htot=Htot,Vtot=Vtot,H0=H0,V0=V0,B0=B0,D0=D0,R0=R0,I0=I0,I_d0=I_d0,I_ac0=I_ac0,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,population=population)

simulation.integr_sci(0,tsim,0.1)

"""
  Simulation Analysis: Plots

"""

# SEIR:
plt.plot(simulation.t,simulation.S, label = 'S',color = 'blue' )
plt.plot(simulation.t,simulation.E, label = 'E',color = 'cyan')
plt.plot(simulation.t,simulation.I, label = 'I',color = 'red')
plt.plot(simulation.t,simulation.R, label = 'R',color = 'green')
plt.plot(simulation.t,simulation.B, label = 'B',color = 'purple')
plt.legend(loc=0)
plt.title('SEIRD')
plt.show()

# Infected:
plt.plot(simulation.t,simulation.Ias, label = 'Ias',color = 'blue' )
plt.plot(simulation.t,simulation.Imi, label = 'Imi',color = 'cyan')
plt.plot(simulation.t,simulation.Ise, label = 'Ise',color = 'red')
plt.plot(simulation.t,simulation.Icr, label = 'Icr',color = 'green')
plt.legend(loc=0)
plt.title('Active Infected')
plt.show()

# Infected:
plt.plot(simulation.t,simulation.Ias_d, label = 'Ias',color = 'blue' )
plt.plot(simulation.t,simulation.Imi_d, label = 'Imi',color = 'cyan')
plt.plot(simulation.t,simulation.Ise_d, label = 'Ise',color = 'red')
plt.plot(simulation.t,simulation.Icr_d, label = 'Icr',color = 'green')
plt.legend(loc=0)
plt.title('Daily Infected')
plt.show()

# Infected:
plt.plot(simulation.t,simulation.Ias_ac, label = 'Ias',color = 'blue' )
plt.plot(simulation.t,simulation.Imi_ac, label = 'Imi',color = 'cyan')
plt.plot(simulation.t,simulation.Ise_ac, label = 'Ise',color = 'red')
plt.plot(simulation.t,simulation.Icr_ac, label = 'Icr',color = 'green')
plt.legend(loc=0)
plt.title('Accumulated Infected')
plt.show()


# Total Population:
plt.hlines(simulation.population,0,simulation.t[-1],label='Total Population') 
plt.plot(simulation.t,simulation.S+simulation.E+simulation.I+simulation.R+simulation.Hse+simulation.Hout+simulation.V+simulation.D+simulation.B, label = 'Simulation Population',color = 'blue')

plt.legend(loc=0)
plt.title('Population conservation test')
plt.show()

# Hospitalized:
plt.plot(simulation.t,simulation.Hse, label = 'Hse',color = 'red' )
plt.plot(simulation.t,simulation.Hout, label = 'Hout',color = 'green')
plt.plot(simulation.t,simulation.V, label = 'V',color = 'purple')
plt.plot(simulation.t,simulation.Hse +simulation.Hout, label = 'Hse+Hout',color = 'black' )

plt.legend(loc=0)
plt.title('Hospitalized')
plt.show()



# Hospitalized + Deaths:
plt.plot(simulation.t,simulation.H_sat, label = 'H_sat',color = 'red' )
plt.plot(simulation.t,simulation.V_sat, label = 'V_sat',color = 'green')

plt.plot(simulation.t,simulation.Hse_D_d, label = 'Hse_D',color = 'blue' )
plt.plot(simulation.t,simulation.V_D_d, label = 'V_D')
plt.plot(simulation.t,simulation.Ise_D_d, label = 'Imi_D',color = 'cyan')
plt.plot(simulation.t,simulation.Icr_D_d, label = 'Ise_D')

plt.legend(loc=0)
plt.title('Hospitalized to Deaths')
plt.show()


"""

Simulation vs data


"""
initdate = datetime(2020,5,15)
currentdate = datetime.now()
currentday = (currentdate - initdate).days

tstate = '13'

# Import Data
RM = ImportData(tstate=tstate,initdate = initdate)
RM.importdata()

# Simulation parameters
tsim = 1000
beta = 0.2
mu = 1.5
k = 0

SeroPrevFactor=1
expinfection=0


# Quarantines 
max_mob = 0.8
rem_mob = 0.5

#alpha = Quarantine(rem_mob,max_mob=max_mob,qp=0,iqt=0,fqt=1000,movfunct = 'once').alpha(t)
alpha = Quarantine(rem_mob).alpha


# Run Simulation
simulation = SEIRHVD(tsim,beta,mu,alpha,k=k,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,InitialConditions = RM)
simulation.integr_sci(0,tsim,0.1)

# ----------------
# Analyse results
# ----------------
days = currentday+20
index = np.searchsorted(simulation.t,days)

# Active Infected
tr_index = np.searchsorted(RM.tr,days)
plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.I[:index],label='Total Active Infected')
plt.scatter(RM.tr,RM.Ir,label='Real Data')
plt.legend(loc=0)
plt.show()

# Daily Infected
plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.I_d[:index],label='Total Daily new Infected')
plt.scatter(RM.I_d_r_tr,RM.I_d_r,label='Real Data')
plt.legend(loc=0)
plt.show()


# Accumulated Infected
plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.I_ac[:index],label='Total Acummulated Infected')
plt.scatter(RM.I_ac_r_tr,RM.I_ac_r,label='Real Data')
plt.legend(loc=0)
plt.show()

# Deaths
plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.B[:index],label='Total Acummulated Deaths')
plt.scatter(RM.Br_tr,RM.Br,label='Real Data')
plt.legend(loc=0)
plt.show()

# UCI/UTI
plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.Hse[:index]+simulation.Hout[:index],label='UCI/UTI Beds')
plt.scatter(RM.sochimi_tr,RM.Hr,label='Real Data')
plt.legend(loc=0)
plt.show()

# VMI
plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.V[:index],label='VMI Usage')
plt.scatter(RM.sochimi_tr,RM.Vr,label='Real Data')
plt.legend(loc=0)
plt.show()



# Data Analyisis

# Infected Error 
# Actives:
idx = np.searchsorted(simulation.t,RM.tr)
simulation.I
 
res = LA.norm(RM.Ir-simulation.I[idx])

Err = np.sum(abs(RM.Ir-simulation.I[idx]))/np.mean(RM.Ir)




tr_index = np.searchsorted(RM.tr,days)


plt.scatter(RM.tr,RM.Ir,label='Real Data')

plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.I_det[:index],label='Infected')
plt.legend(loc=0)
plt.title('Detectable Active Infected')
plt.show()


plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.I_d_det[:index],label='Daily new Infected')
plt.scatter(RM.I_d_r_tr,RM.I_d_r,label='Real Data')
plt.legend(loc=0)
plt.title('Detected Daily infected')
plt.show()

plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.I_ac_det[:index],label='Total Acummulated Infected')
plt.scatter(RM.I_ac_r_tr,RM.I_ac_r,label='Real Data')
plt.legend(loc=0)
plt.show()



plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.Hse[:index]+simulation.Hout[:index],label='UCI/UTI Beds')
plt.scatter(RM.sochimi_tr,RM.Hr,label='Real Data')
plt.scatter(RM.sochimi_tr,RM.Hr_tot,label='Real Data')
plt.legend(loc=0)
plt.show()


plt.xlim(0,days)
plt.plot(simulation.t[:index],simulation.V[:index],label='VMI Usage')
plt.scatter(RM.sochimi_tr,RM.Vr,label='Real Data')
plt.legend(loc=0)
plt.show()





fig, axs = plt.subplots(3, 2)
axs[0,0].set_xlim(0,days)
axs[0,0].scatter(RM.tr,RM.Ir,label='Real Active Data')
axs[0,0].plot(simulation.t[:index],simulation.I_det[:index],label='Infected')
axs[0,0].set_title('Detectable Active Infected')
axs[0,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[0, 1].set_xlim(0,days)
axs[0, 1].plot(simulation.t[:index],simulation.I_d_det[:index],label='Daily new Infected')
axs[0, 1].scatter(RM.I_d_r_tr,RM.I_d_r,label='Daily Real Data')
#axs[0, 1].legend(loc=0)
axs[0, 1].set_title('Detected Daily infected')
axs[0,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')


axs[1,0].set_xlim(0,days)
axs[1,0].plot(simulation.t[:index],simulation.I_ac_det[:index],label='Detected Acummulated Infected')
axs[1,0].scatter(RM.I_ac_r_tr,RM.I_ac_r,label='Accumulated Real Data')
axs[1,0].set_title('Detected Accumulated Infected')
axs[1,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[1,1].set_xlim(0,days)
axs[1,1].plot(simulation.t[:index],simulation.B[:index],label='Total Acummulated Deaths')
axs[1,1].scatter(RM.Br_tr,RM.Br,label='Real Data')
axs[1,1].set_title('Acmulated Deaths')
axs[1,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[2,0].set_xlim(0,days)
axs[2,0].plot(simulation.t[:index],simulation.Hse[:index]+simulation.Hout[:index],label='UCI/UTI Beds')
axs[2,0].scatter(RM.sochimi_tr,RM.Hr,label='Real Data')
axs[2,0].scatter(RM.sochimi_tr,RM.Hr_tot,label='Capacity Data', color = 'red')
axs[2,0].set_title('UCI/UTI Usage')
axs[2,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[2,1].set_xlim(0,days)
axs[2,1].plot(simulation.t[:index],simulation.V[:index],label='VMI Usage')
axs[2,1].scatter(RM.sochimi_tr,RM.Vr,label='Real Data')
axs[2,1].scatter(RM.sochimi_tr,RM.Vr_tot,label='Capacity Data',color = 'red')
axs[2,1].set_title('UCI/UTI Usage')
axs[2,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')




for ax in axs.flat:
    ax.set(xlabel='x-label', ylabel='y-label')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

t_end_idx = np.where(np.array(RM.tr)>=t_end)[0][0]
idx = np.searchsorted(simulation.t,RM.tr[:t_end_idx])
E_I = np.sum(abs(RM.Ir[:t_end_idx]-simulation.I_det[idx]))/(np.mean(RM.Ir[:t_end_idx])*t_end_idx)

t_end_idx = np.where(np.array(RM.I_d_r_tr)>=t_end)[0][0]
idx = np.searchsorted(simulation.t,RM.I_d_r_tr[:t_end_idx])
E_Id = np.sum(abs(RM.I_d_r[:t_end_idx]-simulation.I_d_det[idx]))/(np.mean(RM.I_d_r[:t_end_idx])*t_end_idx)

t_end_idx = np.where(np.array(RM.I_ac_r_tr)>=t_end)[0][0]
idx = np.searchsorted(simulation.t,RM.I_ac_r_tr[:t_end_idx])
E_Iac = np.sum(abs(RM.I_ac_r[:t_end_idx]-simulation.I_ac_det[idx]))/(np.mean(RM.I_ac_r[:t_end_idx])*t_end_idx)

t_end_idx = np.where(np.array(RM.Br_tr)>=t_end)[0][0]
idx = np.searchsorted(simulation.t,RM.Br_tr[:t_end_idx])
E_D = np.sum(abs(RM.Br[:t_end_idx]-simulation.B[idx]))/(np.mean(RM.Br[:t_end_idx])*t_end_idx)

t_end_idx = np.where(np.array(RM.sochimi_tr)>=t_end)[0][0]
idx = np.searchsorted(simulation.t,RM.sochimi_tr[:t_end_idx])
E_UCI = np.sum(abs(RM.Hr[:t_end_idx]-(simulation.Hse[idx]+simulation.Hout[idx])))/(np.mean(RM.Hr[:t_end_idx])*t_end_idx)

t_end_idx = np.where(np.array(RM.sochimi_tr)>=t_end)[0][0]
idx = np.searchsorted(simulation.t,RM.sochimi_tr[:t_end_idx])
E_VMI = np.sum(abs(RM.Vr[:t_end_idx]-simulation.V[idx]))/(np.mean(RM.Vr[:t_end_idx])*t_end_idx)




import math
math.isnan()
E_UCI = np.sum(abs(RM.Hr[:t_end_idx]-(simulation.Hse[idx]+simulation.Hout[idx])))/(np.mean(RM.Hr[:t_end_idx])*t_end_idx)





fig, axs = plt.subplots(3, 2)
axs[0,0].set_xlim(0,days)
axs[0,0].scatter(RM.tr,RM.Ir,label='Real Active Data')
axs[0,0].plot(simulation.t[:index],simulation.I_det[:index],label='Infected')
axs[0,0].set_title('Detectable Active Infected')
axs[0,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[0, 1].set_xlim(0,days)
axs[0, 1].plot(simulation.t[:index],simulation.I_d_det[:index],label='Daily new Infected')
axs[0, 1].scatter(RM.I_d_r_tr,RM.I_d_r,label='Daily Real Data')
#axs[0, 1].legend(loc=0)
axs[0, 1].set_title('Detected Daily infected')
axs[0,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')


axs[1,0].set_xlim(0,days)
axs[1,0].plot(simulation.t[:index],simulation.I_ac_det[:index],label='Detected Acummulated Infected')
axs[1,0].scatter(RM.I_ac_r_tr,RM.I_ac_r,label='Accumulated Real Data')
axs[1,0].set_title('Detected Accumulated Infected')
axs[1,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[1,1].set_xlim(0,days)
axs[1,1].plot(simulation.t[:index],simulation.B[:index],label='Total Acummulated Deaths')
axs[1,1].scatter(RM.Br_tr,RM.Br,label='Real Data')
axs[1,1].set_title('Acmulated Deaths')
axs[1,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[2,0].set_xlim(0,days)
axs[2,0].plot(simulation.t[:index],simulation.Hse[:index]+simulation.Hout[:index],label='UCI/UTI Beds')
axs[2,0].scatter(RM.sochimi_tr,RM.Hr,label='Real Data')
axs[2,0].scatter(RM.sochimi_tr,RM.Hr_tot,label='Capacity Data', color = 'red')
axs[2,0].set_title('UCI/UTI Usage')
axs[2,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[2,1].set_xlim(0,days)
axs[2,1].plot(simulation.t[:index],simulation.V[:index],label='VMI Usage')
axs[2,1].scatter(RM.sochimi_tr,RM.Vr,label='Real Data')
axs[2,1].scatter(RM.sochimi_tr,RM.Vr_tot,label='Capacity Data',color = 'red')
axs[2,1].set_title('UCI/UTI Usage')
axs[2,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')



simulation.pE_Ias = pE_Ias
simulation.pE_Imi = pE_Imi
simulation.pE_Ise = pE_Ise
simulation.pE_Icr = pE_Icr



fig, axs = plt.subplots(3, 2)
axs[0,0].set_xlim(0,days)
axs[0,0].scatter(RM.tr,RM.Ir,label='Real Active Data')
axs[0,0].plot(simulation1.t[:index1],simulation1.I_det[:index1],label='Infected')
axs[0,0].plot(simulation2.t[:index2],simulation2.I_det[:index2],label='Infected')
axs[0,0].plot(simulation3.t[:index3],simulation3.I_det[:index3],label='Infected')
axs[0,0].set_title('Detectable Active Infected')
axs[0,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[0, 1].set_xlim(0,days)
axs[0, 1].plot(simulation1.t[:index1],simulation1.I_d_det[:index1],label='Daily new Infected')
axs[0, 1].plot(simulation2.t[:index2],simulation2.I_d_det[:index2],label='Daily new Infected')
axs[0, 1].plot(simulation3.t[:index3],simulation3.I_d_det[:index3],label='Daily new Infected')
axs[0, 1].scatter(RM.I_d_r_tr,RM.I_d_r,label='Daily Real Data')
#axs[0, 1].legend(loc=0)
axs[0, 1].set_title('Detected Daily infected')
axs[0,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')


axs[1,0].set_xlim(0,days)
axs[1,0].plot(simulation1.t[:index1],simulation1.I_ac_det[:index1],label='Detected Acummulated Infected')
axs[1,0].plot(simulation2.t[:index2],simulation2.I_ac_det[:index2],label='Detected Acummulated Infected')
axs[1,0].plot(simulation3.t[:index3],simulation3.I_ac_det[:index3],label='Detected Acummulated Infected')
axs[1,0].scatter(RM.I_ac_r_tr,RM.I_ac_r,label='Accumulated Real Data')
axs[1,0].set_title('Detected Accumulated Infected')
axs[1,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[1,1].set_xlim(0,days)
axs[1,1].plot(simulation1.t[:index1],simulation1.B[:index1],label='Total Acummulated Deaths')
axs[1,1].plot(simulation2.t[:index2],simulation2.B[:index2],label='Total Acummulated Deaths')
axs[1,1].plot(simulation3.t[:index3],simulation3.B[:index3],label='Total Acummulated Deaths')
axs[1,1].scatter(RM.Br_tr,RM.Br,label='Real Data')
axs[1,1].set_title('Acmulated Deaths')
axs[1,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[2,0].set_xlim(0,days)
axs[2,0].plot(simulation1.t[:index1],simulation1.Hse[:index1]+simulation1.Hout[:index1],label='UCI/UTI Beds')
axs[2,0].plot(simulation2.t[:index2],simulation2.Hse[:index2]+simulation2.Hout[:index2],label='UCI/UTI Beds')
axs[2,0].plot(simulation3.t[:index3],simulation3.Hse[:index3]+simulation3.Hout[:index3],label='UCI/UTI Beds')
axs[2,0].scatter(RM.sochimi_tr,RM.Hr,label='Real Data')
axs[2,0].scatter(RM.sochimi_tr,RM.Hr_tot,label='Capacity Data', color = 'red')
axs[2,0].set_title('UCI/UTI Usage')
axs[2,0].axvline(x = t_end, linestyle = 'dotted',color = 'grey')

axs[2,1].set_xlim(0,days)
axs[2,1].plot(simulation1.t[:index1],simulation1.V[:index1],label='VMI Usage')
axs[2,1].plot(simulation2.t[:index2],simulation2.V[:index2],label='VMI Usage')
axs[2,1].plot(simulation3.t[:index3],simulation3.V[:index3],label='VMI Usage')
axs[2,1].scatter(RM.sochimi_tr,RM.Vr,label='Real Data')
axs[2,1].scatter(RM.sochimi_tr,RM.Vr_tot,label='Capacity Data',color = 'red')
axs[2,1].set_title('UCI/UTI Usage')
axs[2,1].axvline(x = t_end, linestyle = 'dotted',color = 'grey')
