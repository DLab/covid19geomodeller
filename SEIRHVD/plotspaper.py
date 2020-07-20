#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -------------------- #
#                      #
#     SEIRHDV Paper    #
#                      #
# -------------------- #

from SEIRHVD_local import SEIRHVD_local
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import multiprocessing  
from joblib import Parallel, delayed

tstate = ''
initdate = datetime(2020,5,15)

# ------------------- #
#        Plot 1       #
# ------------------- # 

# Camas H totales vs infectados severos vs HcrtoD




# Parametros del modelo
beta = 0.117 # Tasa de contagio
mu = 0.6 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos

tsim = 1000 # Tiempo de simulacion



# Simular
# Activos Iniciales
I_act0 = 100
I_as_prop = 1
I_mi_prop = 0
I_se_prop = 0
I_cr_prop = 0

# Muertos iniciales
dead0 = 0

population = 1000000
# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0
# UCI/UTI capacity per 1000 persons
nm = int(population/100000)



# Hospital capacity
step = 6
Htot_max = 30
Htot = list(range(0*nm,Htot_max*nm+step*nm,step*nm))
Htot_per1000 = list(range(0,30+step,step))

# VMI Capacity
Vtot = 10*nm

#Movilty
# From 0 to 1 in steps of:
step = 0.25
alpha = list(np.arange(0,1+step,step))

Nalpha = len(alpha)
NHtot = len(Htot)


# Simulation

sims = []
for i in Htot:    
    aux = []
    for j in alpha:
        print("alpha")
        # Creación del objeto de simulación 
        simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
        quarantines = [[tsim, 0.85, j, 0.0, 0.0, tsim, 0.0]]
        simulation.inputarray = np.array(quarantines)
        simulation.addquarantine()
        simulation.initialvalues(I_act0,dead0,population,H0,V0,i,Vtot,R=0,D=0,H_cr = 0)
        simulation.simulate(v=3)
        if i > 0:
            aux.append(simulation)
        else:
            aux.append(simulation)
    sims.append(aux)

# Plots Generation:
maxval = max([max([max(sims[i][j].D[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = 600
fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(NHtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_crD_d[0],label="H_cr to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].VD_d[0],label="V to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_crD_d[0],label="I_cr to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_seD_d[0],label="I_se to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].D[0],label="Deaths")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])
fig.suptitle('Axes values are scaled individually by default')
#fig.tight_layout()
lines, labels = fig.axes[-1].get_legend_handles_labels()
    
#fig.legend(lines, labels, loc = 'upper center')
fig.legend(lines, labels,loc = 'best')
fig.show()

# Cambiar X to D Diarios 

# ------------------- #
#        Plot 2       #
# ------------------- # 

# Ventiladores vs infectados criticos vs IcrtoD vs VtoD 



# Parametros del modelo
beta = 0.117 # Tasa de contagio
mu = 0.6 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos

tsim = 1000 # Tiempo de simulacion




# Activos Iniciales
I_act0 = 100
I_as_prop = 1
I_mi_prop = 0
I_se_prop = 0
I_cr_prop = 0

# Muertos iniciales
dead0 = 0

population = 1000000
# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0
# UCI/UTI capacity per 100.000 persons
nm = int(population/100000)
step = 6*nm

Htot = 30#list(range(0*nm,30*nm+step,step))

# VMI Capacity
Vtot_max = 20 #per 100.000 persons
Vtot = list(range(0*nm,Vtot_max*nm+step,step))

#Movilty
# From 0 to 1 in steps of:
step = 0.25
alpha = list(np.arange(0,1+step,step))

Nalpha = len(alpha)
NVtot = len(Vtot)


# Simulation

sims = []
for i in Vtot:    
    aux = []
    for j in alpha:
        print("alpha")
        # Creación del objeto de simulación 
        simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)
        quarantines = [[tsim, 0.85, j, 0.0, 0.0, tsim, 0.0]]
        simulation.inputarray = np.array(quarantines)
        simulation.addquarantine()
        simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,i,R=0,D=0,H_cr = 0,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
        simulation.simulate(v=3)
        if i > 0:
            aux.append(simulation)
        else:
            aux.append(simulation)
    sims.append(aux)

# Plots Generation:

fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_crD[0],label="H_cr to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].VD[0],label="V to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_crD[0],label="I_cr to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_seD[0],label="I_se to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].B[0],label="Deaths")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
fig.suptitle('Axes values are scaled individually by default')
#fig.tight_layout()
lines, labels = fig.axes[-1].get_legend_handles_labels()
    
#fig.legend(lines, labels, loc = 'upper center')
fig.legend(lines, labels,loc = 'best')
fig.show()





# ------------------- #
#        Plot 3       #
# ------------------- # 

# Contourplot de SHFR=Muertos acumulados/((Ise+Icr) acumulados) considerando movilidad (alpha) en el eje X, 
# y numero de camas en el eje Y
tstate = ''
initdate = datetime(2020,5,15)


# Parametros del modelo
beta = 0.12 # Tasa de contagio
mu = 0.6 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos

tsim = 1000 # Tiempo de simulacion


# Simular
# Activos Iniciales
I_act0 = 100
I_as_prop = 1
I_mi_prop = 0
I_se_prop = 0
I_cr_prop = 0



# Muertos iniciales
dead0 = 0

population = 100000
# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0
# UCI/UTI capacity per 100000 persons
Htot_max = 50
nm = int(population/100000)
step = 2
Htot = list(range(0*nm,Htot_max*nm+step*nm,step*nm))
Htot_per100M = list(range(0,Htot_max+step,step))
# VMI Capacity
Vtot = [i/2 for i in Htot]

# Movility
step = 0.05
alpha = list(np.arange(0.05,1+step,step))


# Simulation 
num_cores = multiprocessing.cpu_count() 
def SHFRsimulate(Htot,Vtot,alpha): 
    simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop) 
    quarantines = [[tsim, 0.85, alpha, 0.0, 0.0, tsim, 0.0]] 
    simulation.inputarray = np.array(quarantines) 
    simulation.addquarantine() 
    simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
    simulation.simulate(v=3)  
    return simulation #, simulation.SHFR[0]  
 
# Run  Simulation 
SHFR = np.zeros((len(Htot),len(alpha)))   
sims = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(SHFRsimulate)(Htot[i],Vtot[i],alpha[j]) for j in range(len(alpha)))     
    sims.append(aux)                 
 
for i in range(len(Htot)):     
    for j in range(len(alpha)): 
        SHFR[i][j] = sims[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,Htot_per100M,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 


# Non parallel Simulation
SHFR = np.zeros((len(Htot),len(alpha)))  
sims = []
for i in range(len(Htot)):
    aux = []
    for j in range(len(alpha)):
        # Creación del objeto de simulación 
        simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)
        quarantines = [[tsim, 0.85, alpha[j], 0.0, 0.0, tsim, 0.0]]
        simulation.inputarray = np.array(quarantines)
        simulation.addquarantine()
        simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot[i],Vtot,R=0,D=0,H_cr = 0,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
        simulation.simulate(v=3)
        SHFR[i,j] = simulation.SHFR[0]
        aux.append(simulation)
    sims.append(aux)        





# SHFR QA:
tstate = ''
initdate = datetime(2020,5,15)


# Parametros del modelo
beta = 0.2 # Tasa de contagio
mu = 0.6 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos

tsim = 1000 # Tiempo de simulacion


# Simular
# Activos Iniciales
I_act0 = 100
I_as_prop = 1
I_mi_prop = 0
I_se_prop = 0
I_cr_prop = 0

# Muertos iniciales
dead0 = 0

population = 100000
# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0
# UCI/UTI capacity per 100000 persons
Htot = 50

# VMI Capacity
Vtot = 25
# Movility
alpha = 0.8

simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)
quarantines = [[tsim, 0.85, alpha, 0.0, 0.0, tsim, 0.0]]
simulation.inputarray = np.array(quarantines)
simulation.addquarantine()
simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
simulation.simulate(v=3)

# Análisis de Resultados