#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# -------------------- #
#                      #
#     SEIRHDV Paper    #
#                      #
# -------------------- #

Generación de Plots Finales
- beta: 0.2
- init: infectados 100 todos asintomaticos
- mu: 0 
- 

- Movilidad basal 0.3
- Movilidad superior: 0.65,0.70, 0.75
- Cuarentena total, cuarentena dinamica 2 semanas
- mass action, saturated kinetics con SKF 30

Grillas
- Camas
- Ventiladores
- SHFR(t)

Contour:
- SHFR
- MassAction vs SK

"""
import sys
from pathlib import Path
sys.path.insert(1, '../src/')
sys.path.insert(1, 'src/')
from SEIRHVD_local import SEIRHVD_local
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import multiprocessing  
from joblib import Parallel, delayed
from numpy import linalg as LA
import os

tstate = ''
initdate = datetime(2020,5,15)

beep = lambda x: os.system("echo -n '\a';sleep 0.2;" * x)

num_cores = multiprocessing.cpu_count() 
def ParallelSimulation(Htot,Vtot,max_mov = 0.85,rem_mov = 0.65,k=0, qp = 0, qt = 0,iqt = 0,fqt = 500): 
    simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop, k = k) 
    quarantines = [[tsim, max_mov, rem_mov, qp, iqt, fqt, qt]] 
    simulation.inputarray = np.array(quarantines) 
    simulation.addquarantine() 
    simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
    simulation.simulate(v=3)  
    return simulation #, simulation.SHFR[0]  


# ------------------------------------- #
#                                       #
#     Parametros Epi Generales Paper    #
#                                       #
# ------------------------------------- #

# Parametros del modelo
beta = 0.2 # Tasa de contagio
mu = 0 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos

# Activos Iniciales
I_act0 = 100
I_as_prop = 1
I_mi_prop = 0
I_se_prop = 0
I_cr_prop = 0

# Muertos iniciales
dead0 = 0
population = 1000000
nm = int(population/100000)

# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0

# Tiempo de simulacion
tsim = 1000 

# Rem Mov:
rem_mov = 0.3

# Max Mov:
max_mov = [0.65,0.7,0.75]



""" 
 # ------------------------------------------------- #
 #                                                   #
 #     Análisis para total de Camas Disponibles      #
 #                                                   #
 # ------------------------------------------------- #

    Datos: Fallecimientos relacionados a escaséz de camas
    ejes: Camas H totales vs infectados severos vs 

"""

# Hospital capacity
step = 8
Htot_max = 40
Htot = list(range(0*nm,Htot_max*nm+step*nm,step*nm))
Htot_per1000 = list(range(0,30+step,step))
NHtot = len(Htot)

# VMI Capacity
Vtot = 10*nm

#Movilty
# From 0 to 1 in steps of:
max_mov = [0.3,0.65,0.75,0.85]
rem_mov = 0.3
Nalpha = len(max_mov)


# ---------------------------------------------------- #
#         Cuarentena Total   -   Mass Action           #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim


# Simulation

sims1 = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims1.append(aux)


# Plots Generation:
maxval = max([max([max(sims1[i][j].I_se[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = 300


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        #axs[i, j].plot(sims1[i][j].t[0],sims1[i][j].H_bed[0]+sims1[i][j].I_se[0]+sims1[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims1[i][j].t[0],sims1[i][j].H_bed[0],label="Hospitalized")
        axs[i, j].plot(sims1[i][j].t[0],sims1[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims1[i][j].t[0],sims1[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot per 100k : "+str(int(Htot[i]/10))+" | Mov: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Bed Analysis - Total Quarantine')
fig.show()

# Analysis
I_se_peak1 = np.zeros((NHtot,Nalpha)) 
I_seD_peak1 = np.zeros((NHtot,Nalpha)) 
peak_day1 = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):
    for j in range(Nalpha):
        I_se_peak1[i][j] = max(sims1[i][j].I_se[0])
        I_seD_peak1[i][j] = max(sims1[i][j].I_seD_d[0])
        peak_day1[i][j] = sims1[i][j].peak_t[0]



# ---------------------------------------------------- #
#        Cuarentena Total   - Saturated kinetics       #
# ---------------------------------------------------- # 

# Mass Action
k = 30
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim


# Simulation

sims2 = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims2.append(aux)
    

# Plots Generation:
maxval = max([max([max(sims2[i][j].I_se[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = 300


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        #axs[i, j].plot(sims2[i][j].t[0],sims2[i][j].H_bed[0]+sims2[i][j].I_se[0]+sims2[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims2[i][j].t[0],sims2[i][j].H_bed[0],label="Hospitalized")
        axs[i, j].plot(sims2[i][j].t[0],sims2[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims2[i][j].t[0],sims2[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot per 100k : "+str(int(Htot[i]/10))+" | Mov: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Bed Analysis - Total Quarantine - SKF=30')
fig.show()

# Analysis
I_se_peak2 = np.zeros((NHtot,Nalpha)) 
I_seD_peak2 = np.zeros((NHtot,Nalpha)) 
peak_day2 = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):
    for j in range(Nalpha):
        I_se_peak2[i][j] = max(sims2[i][j].I_se[0])
        I_seD_peak2[i][j] = max(sims2[i][j].I_seD_d[0])
        peak_day2[i][j] = sims2[i][j].peak_t[0]



# --------------------------------------------------------------- #
#        Cuarentena Dinamica 14 dias -   Mass Action              #
# --------------------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim


# Simulation

sims3 = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,max_mov = j, rem_mov = rem_mov,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims3.append(aux)
    

# Plots Generation:
maxval = max([max([max(sims3[i][j].I_se[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = 500


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        #axs[i, j].plot(sims3[i][j].t[0],sims3[i][j].H_bed[0]+sims3[i][j].I_se[0]+sims3[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims3[i][j].t[0],sims3[i][j].H_bed[0],label="Hospitalized")
        axs[i, j].plot(sims3[i][j].t[0],sims3[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims3[i][j].t[0],sims3[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot per 100k : "+str(int(Htot[i]/10))+" | Mov: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Bed Analysis - Total Quarantine')
fig.show()

# Analysis
I_se_peak3 = np.zeros((NHtot,Nalpha)) 
I_seD_peak3 = np.zeros((NHtot,Nalpha)) 
peak_day3 = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):
    for j in range(Nalpha):
        I_se_peak3[i][j] = max(sims3[i][j].I_se[0])
        I_seD_peak3[i][j] = max(sims3[i][j].I_seD_d[0])
        peak_day3[i][j] = sims3[i][j].peak_t[0]





# --------------------------------------------------------------- #
#        Cuarentena Dinamica 14 dias - Saturated kinetics         #
# --------------------------------------------------------------- # 

# Saturated kinetics Factor
k = 30
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim


# Simulation

sims4 = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,max_mov = j, rem_mov = rem_mov,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims4.append(aux)
    




# Plots Generation:
maxval = max([max([max(sims4[i][j].I_se[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = 1000


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        #axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].H_bed[0]+sims4[i][j].I_se[0]+sims4[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].H_bed[0],label="Hospitalized")
        axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot per 100k : "+str(int(Htot[i]/10))+" | Mov: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Bed Analysis - Total Quarantine - SKF = 30')
fig.show()

# Analysis
I_se_peak4 = np.zeros((NHtot,Nalpha)) 
I_seD_peak4 = np.zeros((NHtot,Nalpha)) 
peak_day4 = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):
    for j in range(Nalpha):
        I_se_peak4[i][j] = max(sims4[i][j].I_se[0])
        I_seD_peak4[i][j] = max(sims4[i][j].I_seD_d[0])
        peak_day4[i][j] = sims4[i][j].peak_t[0]







""" 
 # ------------------------------------------------- #
 #                                                   #
 #      Análisis para total de VMI Disponibles       #
 #                                                   #
 # ------------------------------------------------- #

    Datos: Fallecimientos relacionados a escaséz de VMI
    ejes: VMI totales vs infectados severos 

"""
Htot = 30#list(range(0*nm,30*nm+step,step))

# VMI Capacity
Vtot_max = 20 #per 100.000 persons
step = 4*nm
Vtot = list(range(0*nm,Vtot_max*nm+step,step))
NVtot = len(Vtot)

#Mobilty
# From 0 to 1 in steps of:
max_mov = [0.3,0.65,0.75,0.85]
rem_mov = 0.3
Nalpha = len(max_mov)


# ---------------------------------------------------- #
#         Cuarentena Total   -   Mass Action           #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Simulation
sims5 = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims5.append(aux)

# Plots Generation:

maxval = max([max([max(sims5[i][j].I_cr[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = 350

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('VMI Analysis - Total Quarantine')
fig.show()


# Analysis
I_cr_peak5 = np.zeros((NHtot,Nalpha)) 
I_crD_peak5 = np.zeros((NHtot,Nalpha)) 
peak_day5 = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):
    for j in range(Nalpha):
        I_cr_peak5[i][j] = max(sims5[i][j].I_cr[0])
        I_crD_peak5[i][j] = max(sims5[i][j].I_crD_d[0])
        peak_day5[i][j] = sims5[i][j].peak_t[0]


# ---------------------------------------------------- #
#       Cuarentena Total   -   Saturated kinetics      #
# ---------------------------------------------------- # 

# Mass Action
k = 30
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Simulation
sims6 = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims6.append(aux)

# Plots Generation:
maxval = max([max([max(sims6[i][j].I_cr[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('VMI Analysis - Total Quarantine - SKF = 30')
fig.show()

# Analysis
I_cr_peak6 = np.zeros((NHtot,Nalpha)) 
I_crD_peak6 = np.zeros((NHtot,Nalpha)) 
peak_day6 = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):
    for j in range(Nalpha):
        I_cr_peak6[i][j] = max(sims6[i][j].I_cr[0])
        I_crD_peak6[i][j] = max(sims6[i][j].I_crD_d[0])
        peak_day6[i][j] = sims6[i][j].peak_t[0]


# ---------------------------------------------------- #
#    Dynamic Quarantine 14 days  -   Mass Action       #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Simulation
sims7 = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,max_mov = j, rem_mov = rem_mov,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims7.append(aux)

# Plots Generation:

maxval = max([max([max(sims7[i][j].I_cr[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = 500

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('VMI Analysis - Dynamic Quarantine  14 days')
fig.show()


# Analysis
I_cr_peak7 = np.zeros((NHtot,Nalpha)) 
I_crD_peak7 = np.zeros((NHtot,Nalpha)) 
peak_day7 = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):
    for j in range(Nalpha):
        I_cr_peak7[i][j] = max(sims7[i][j].I_cr[0])
        I_crD_peak7[i][j] = max(sims7[i][j].I_crD_d[0])
        peak_day7[i][j] = sims7[i][j].peak_t[0]






# ---------------------------------------------------------------- #
#        Dynamic Quarantine 14 days  -   Saturated kinetics        #
# ---------------------------------------------------------------- # 

# Mass Action
k = 30
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Simulation
sims8 = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,max_mov = j, rem_mov = rem_mov,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims8.append(aux)

# Plots Generation:

maxval = max([max([max(sims8[i][j].I_cr[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims8[i][j].t[0],sims8[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims8[i][j].t[0],sims8[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims8[i][j].t[0],sims8[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims8[i][j].t[0],sims8[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims8[i][j].t[0],sims8[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('VMI Analysis - Dynamic Quarantine  14 days - SKF = 30')
fig.show()

# Analysis
I_cr_peak6 = np.zeros((NHtot,Nalpha)) 
I_crD_peak6 = np.zeros((NHtot,Nalpha)) 
peak_day6 = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):
    for j in range(Nalpha):
        I_cr_peak6[i][j] = max(sims6[i][j].I_cr[0])
        I_crD_peak6[i][j] = max(sims6[i][j].I_crD_d[0])
        peak_day6[i][j] = sims6[i][j].peak_t[0]





""" 
 # ---------------------------------------------------------------------- #
 #                                                                        #
 #                           Análisis SHFR(t)                             #
 #                                                                        #
 # ---------------------------------------------------------------------- #

 """

# Hospital capacity
step = 8
Htot_max = 40
Htot = list(range(0*nm,Htot_max*nm+step*nm,step*nm))
Htot_per1000 = list(range(0,30+step,step))
NHtot = len(Htot)


# VMI Capacity
Vtot = [i/2 for i in Htot]



#Mobilty
# From 0 to 1 in steps of:
max_mov = [0.3,0.65,0.75,0.85]
rem_mov = 0.3
Nalpha = len(max_mov)

# ---------------------------------------------------------------------------- #
#                 SHFR(t) -   Total Quarantine  - SKF = 0                      #
# ---------------------------------------------------------------------------- # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Run  Simulation 
sims9 = []
for i in range(NHtot):
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims9.append(aux)

beep(10)

xlim = 400
ylim = 1

fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(len(Htot)):
    for j in range(Nalpha):
        axs[i, j].plot(sims9[i][j].t[0],sims9[i][j].SHFR_d[0],label="SHFR")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,ylim*1.05])
        axs[i, j].set_xlim([0,xlim])            

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('SHFR(t) - Total Quarantine  - SFK = 0')
fig.show()


# ---------------------------------------------------------------------------- #
#                 SHFR(t) -   Total Quarantine  - SKF = 30                      #
# ---------------------------------------------------------------------------- # 
 
 # Mass Action
k = 30
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Run  Simulation 
sims10 = []
for i in range(NHtot):
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims10.append(aux)

beep(10)

xlim = 400
ylim = 1


fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(len(Htot)):
    for j in range(Nalpha):
        axs[i, j].plot(sims10[i][j].t[0],sims10[i][j].SHFR_d[0],label="SHFR")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,ylim*1.05])
        axs[i, j].set_xlim([0,xlim])             

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('SHFR(t) - Total Quarantine  - SFK = 30')
fig.show()



# ------------------------------------------------------------------------ #
#              SHFR(t) Dynamic Quarantine 14 days - SKF = 0                #
# ------------------------------------------------------------------------ # 


# Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
sims11 = []
for i in range(NHtot):
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],max_mov = j, rem_mov = rem_mov,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims11.append(aux)

beep(10)

xlim = 400
ylim = 1


fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(len(Htot)):
    for j in range(Nalpha):
        axs[i, j].plot(sims11[i][j].t[0],sims11[i][j].SHFR_d[0],label="SHFR")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,ylim*1.05])
        axs[i, j].set_xlim([0,xlim])             

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('SHFR(t) -  Dynamic Quarantine  - SFK = 0')
fig.show()



# ------------------------------------------------------------------------ #
#              SHFR(t) Dynamic Quarantine 14 days - SKF = 30                #
# ------------------------------------------------------------------------ # 


# Mass Action
k = 30
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
sims12 = []
for i in range(NHtot):
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],max_mov = j, rem_mov = rem_mov,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims12.append(aux)

beep(10)

xlim = 400
ylim = 1


fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(len(Htot)):
    for j in range(Nalpha):
        axs[i, j].plot(sims12[i][j].t[0],sims12[i][j].SHFR_d[0],label="SHFR")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(max_mov[j]))
        axs[i, j].set_ylim([0,ylim*1.05])
        axs[i, j].set_xlim([0,xlim])             

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('SHFR(t) - Dynamic Quarantine  - SFK = 30')
fig.show()






""" 
 # ---------------------------------------------------------------------- #
 #                                                                        #
 #                Análisis SHFR final simulacion                          #
 #                                                                        #
 # ---------------------------------------------------------------------- #

    Datos: SHFR
    ejes: VMI totales vs infectados severos 


    Contourplot de SHFR=Muertos acumulados/((Ise+Icr) acumulados) considerando movilidad (alpha) en el eje X, 
     y numero de camas en el eje Y

"""

tsim = 2000
# Hospital capacity
step = 5
Htot_max = 50
Htot = list(range(0*nm,Htot_max*nm+step*nm,step*nm))
Htot_per100k = [int(Htot[i]/nm) for i in range(len(Htot))]
NHtot = len(Htot)


# VMI Capacity
Vtot = [i/2 for i in Htot]


#Mobilty
# From 0 to 1 in steps of:
max_mov = [0.3,0.65,0.75,0.85]
step =0.05
max_mov= list(np.arange(0.05,1+step,step))
rem_mov = 0.3
Nalpha = len(max_mov)


# ---------------------------------------------------------------------------- #
#                 SHFR -   Total Quarantine  - SKF = 0                      #
# ---------------------------------------------------------------------------- # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Run  Simulation 
sims13 = []
for i in range(NHtot):
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims13.append(aux)
  
beep(10)      

SHFR = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):     
    for j in range(Nalpha): 
        SHFR[i][j] = sims13[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(max_mov,Htot_per100k,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR - Total Quarantine - SKF = 0')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 

# ---------------------------------------------------------------------------- #
#                 SHFR -   Total Quarantine  - SKF = 30                      #
# ---------------------------------------------------------------------------- # 
 
 # Mass Action
k = 30
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Run  Simulation 
sims14 = []
for i in range(NHtot):
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims14.append(aux)
  
beep(10)      

SHFR = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):     
    for j in range(Nalpha): 
        SHFR[i][j] = sims14[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(max_mov,Htot_per100k,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR - Total Quarantine - SKF = 30 ')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 





# ---------------------------------------------------------------------------- #
#                 SHFR -   Dynamic Quarantine  14 days - SKF = 0               #
# ---------------------------------------------------------------------------- # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
sims15 = []
for i in range(NHtot):
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims15.append(aux)
  
beep(10)      

SHFR = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):     
    for j in range(Nalpha): 
        SHFR[i][j] = sims15[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(max_mov,Htot_per100k,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR - Dynamic Quarantine 14 days - SKF = 0')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 





# ---------------------------------------------------------------------------- #
#                 SHFR -   Dynamic Quarantine  14 days - SKF = 30               #
# ---------------------------------------------------------------------------- # 
 
 # Mass Action
k = 30
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
sims16 = []
for i in range(NHtot):
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],max_mov = j, rem_mov = j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in max_mov)     
    sims16.append(aux)
  
beep(10)      

SHFR = np.zeros((NHtot,Nalpha)) 
for i in range(NHtot):     
    for j in range(Nalpha): 
        SHFR[i][j] = sims16[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(max_mov,Htot_per100k,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR - Dynamic Quarantine 14 days - SKF = 30')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 










""" 
 # ---------------------------------------------------------------------- #
 #                                                                        #
 #               Análisis Mass Action vs Saturated Kynetics               #
 #                                                                        #
 # ---------------------------------------------------------------------- #

# Comparacion dinamica Mass Action con Saturated Kinetics
# Grafico se calcula para un valor T de tiempo. Debe ser contourplor con ejes x= alpha, y=kappa,
# y color dado por funcion distancia entre Peak de Infectados para los valores de kappa respecto al de kappa =0


"""
tsim = 2000
# UCI/UTI capacity
Htot = 50*nm 
# VMI Capacity
Vtot = Htot/2

# Movility
step = 0.05
max_mov = list(np.arange(0.3,0.85+step,step))

# Saturation Kinetics Factor
k = [0,5,10,20,30,40]#list(np.arange(0,kmax+step,step))


# ------------------------------------------------- #
#        Plot 1: MAD vs SKD -Total Quarantine       #
# ------------------------------------------------- # 

qt = 0
qp = 0
iqt = 0
fqt = tsim
sims17 = [] 

# Run  Simulation 
sims17 = [] 
for i in k: 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,Vtot,max_mov=j,rem_mov = j,k=i,qp = qp, qt = qt, iqt = iqt, fqt =fqt) for j in max_mov)
    sims17.append(aux)

beep(10)    

# Peak size proportion
peak = []
for i in range(len(k)):
    aux = []
    for j in range(len(max_mov)):
        aux.append(sims17[i][j].peak[0]/sims17[0][j].peak[0])
    peak.append(aux)
        

# Contour Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(max_mov,k,peak) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('MAD vs SKD - Total Quarantine')
ax.set_xlabel('Mobility')
ax.set_ylabel('Saturation Dynamics Factor ')
plt.show() 


# Grid plot
fig, axs = plt.subplots(len(k), len(max_mov))
for i in range(len(k)):
    for j in range(len(max_mov)):
        axs[i, j].plot(sims17[i][j].t[0],sims17[i][j].I[0],label="Infected")  
        axs[i, j].set_title("K: "+str(k[i])+" | Alpha: "+str(round(max_mov[j],2)))
fig.suptitle('MAD vs SKD  | Htot = 50  | Total Quarantine')
#fig.tight_layout()
lines, labels = fig.axes[-1].get_legend_handles_labels()
    
#fig.legend(lines, labels, loc = 'upper center')
fig.legend(lines, labels,loc = 'best')
fig.show()



# --------------------------------------------------------- #
#        Plot 2: MAD vs SKD -Dynamic Quarantine 14days      #
# --------------------------------------------------------- # 

qt = 1
qp = 14
iqt = 0
fqt = tsim
sims18 = [] 

# Run  Simulation 
sims18 = [] 
for i in k: 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,Vtot,max_mov=j,rem_mov = rem_mov,k=i,qp = qp, qt = qt, iqt = iqt, fqt =fqt) for j in max_mov)
    sims18.append(aux)

beep(10)    

# Peak size proportion
peak = []
for i in range(len(k)):
    aux = []
    for j in range(len(max_mov)):
        aux.append(sims18[i][j].peak[0]/sims18[0][j].peak[0])
    peak.append(aux)
        

# Contour Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(max_mov,k,peak) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Peak Size Proportion - Dynamic Quarantine 14 Days')
ax.set_xlabel('Mobility')
ax.set_ylabel('Saturation Dynamics Factor')
plt.show() 


# Grid plot
fig, axs = plt.subplots(len(k), len(max_mov))
for i in range(len(k)):
    for j in range(len(max_mov)):
        axs[i, j].plot(sims18[i][j].t[0],sims18[i][j].I[0],label="Infected")  
        axs[i, j].set_title("K: "+str(k[i])+" | Alpha: "+str(max_mov[j]))
fig.suptitle('MAD vs SKD  | Htot = 50  | Dynamic Quarantine 14 Days')
#fig.tight_layout()
lines, labels = fig.axes[-1].get_legend_handles_labels()
    
#fig.legend(lines, labels, loc = 'upper center')
fig.legend(lines, labels,loc = 'best')
fig.show()

