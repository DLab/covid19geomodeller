#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -------------------- #
#                      #
#     SEIRHDV Paper    #
#                      #
# -------------------- #
import sys
from pathlib import Path
sys.path.insert(1, '..src/SEIRHVD/')
sys.path.insert(1, 'src/SEIRHVD/')
from SEIRHVD_local import SEIRHVD_local
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import multiprocessing  
from joblib import Parallel, delayed
from numpy import linalg as LA

tstate = ''
initdate = datetime(2020,5,15)


num_cores = multiprocessing.cpu_count() 
def ParallelSimulation(Htot,Vtot,alpha,k=0, qp = 0, qt = 0,iqt = 0,fqt = 500): 
    simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop, k = k) 
    quarantines = [[tsim, 0.85, alpha, qp, iqt, fqt, qt]] 
    simulation.inputarray = np.array(quarantines) 
    simulation.addquarantine() 
    simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
    simulation.simulate(v=3)  
    return simulation #, simulation.SHFR[0]  




""" 
 # ------------------------------------------------- #
 #                                                   #
 #     Análisis para total de Camas Disponibles      #
 #                                                   #
 # ------------------------------------------------- #

    Datos: Fallecimientos relacionados a escaséz de camas
    ejes: Camas H totales vs infectados severos vs 

"""

# Parametros del modelo
beta = 0.2 # Tasa de contagio
mu = 0.6 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos


# Simular
tsim = 1000 # Tiempo de simulacion
population = 1000000
nm = int(population/100000)

# Activos Iniciales
I_act0 = 100
I_as_prop = 1
I_mi_prop = 0
I_se_prop = 0
I_cr_prop = 0

# Muertos iniciales
dead0 = 0

# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0

# Hospital capacity
step = 6
Htot_max = 30
Htot = list(range(0*nm,Htot_max*nm+step*nm,step*nm))
Htot_per1000 = list(range(0,30+step,step))

# VMI Capacity
Vtot = 10*nm

#Movilty
# From 0 to 1 in steps of:
step = 0.05
alpha = list(np.arange(0.5,0.7+step,step))

Nalpha = len(alpha)
NHtot = len(Htot)


# -------------------------------------------------- #
#        Plot 1:  Cuarentena total permanente        #
# -------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Simulation

sims = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)


# Plots Generation:
maxval = max([max([max(sims[i][j].I_se[0]+sims[i][j].H_bed[0]+sims[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = 500#tsim


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0]+sims[i][j].I_se[0]+sims[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0],label="Hospitalized")
        #axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('MassAction - Total Quarantine')
fig.show()

# ---------------------------------------------------- #
#        Plot 2:  Cuarentena dinamica 2 semanas        #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Simulation

sims = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)


# Plots Generation:
maxval = max([max([max(sims[i][j].I_se[0]+sims[i][j].H_bed[0]+sims[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = tsim


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0]+sims[i][j].I_se[0]+sims[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0],label="Hospitalized")
        #axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('MassAction - Dynamic Quarantine (2W)')
fig.show()

# ------------------------------------------------- #
#        Plot 3:  Cuarentena Total SKF 30           #
# ------------------------------------------------- # 

# Mass Action
k = 30
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Simulation

sims = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)


# Plots Generation:
maxval = max([max([max(sims[i][j].I_se[0]+sims[i][j].H_bed[0]+sims[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = tsim


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0]+sims[i][j].I_se[0]+sims[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0],label="Hospitalized")
        #axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Saturated Kynetics -Total Quarantine')
fig.show()



# ----------------------------------------------------------- #
#        Plot 4:  Cuarentena dinamica 2 semanas SKF 30        # 
# ----------------------------------------------------------- # 

# Mass Action
k = 30
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Simulation

sims = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)


# Plots Generation:
maxval = max([max([max(sims[i][j].I_se[0]+sims[i][j].H_bed[0]+sims[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = tsim


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0]+sims[i][j].I_se[0]+sims[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0],label="Hospitalized")
        #axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Saturated Kynetics - Dynamic Quarantine (2W)')
fig.show()






""" 
 # ------------------------------------------------- #
 #                                                   #
 #      Análisis para total de VMI Disponibles       #
 #                                                   #
 # ------------------------------------------------- #

    Datos: Fallecimientos relacionados a escaséz de VMI
    ejes: VMI totales vs infectados severos 

"""

# Parametros del modelo
beta = 0.2 # Tasa de contagio
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
nm = int(population/100000)

# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0


Htot = 30#list(range(0*nm,30*nm+step,step))

# VMI Capacity
Vtot_max = 20 #per 100.000 persons
step = 6*nm
Vtot = list(range(0*nm,Vtot_max*nm+step,step))

#Mobilty

# From 0 to 1 in steps of:
step = 0.25
alpha = list(np.arange(0.5,0.7+step,step))
Nalpha = len(alpha)
NVtot = len(Vtot)


# -------------------------------------------------- #
#        Plot 1:  Cuarentena total permanente        #
# -------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Simulation

sims = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)

# Plots Generation:

maxval = max([max([max(sims[i][j].I_cr[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('MassAction - Total Quarantine')
fig.show()


# ---------------------------------------------------- #
#        Plot 2:  Cuarentena dinamica 2 semanas        #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Simulation

sims = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)

# Plots Generation:

maxval = max([max([max(sims[i][j].I_se[0]+sims[i][j].H_bed[0]+sims[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(alpha[j]))
        #axs[i, j].set_ylim([0,maxval*1.05])
        #axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('MassAction - Total Quarantine')
fig.show()

# ------------------------------------------------- #
#        Plot 3:  Cuarentena Total SKF 30           #
# ------------------------------------------------- # 

# Mass Action
k = 30
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Simulation

sims = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)

# Plots Generation:

maxval = max([max([max(sims[i][j].I_se[0]+sims[i][j].H_bed[0]+sims[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(alpha[j]))
        #axs[i, j].set_ylim([0,maxval*1.05])
        #axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('MassAction - Total Quarantine')
fig.show()


# ----------------------------------------------------------- #
#        Plot 4:  Cuarentena dinamica 2 semanas SKF 30        # 
# ----------------------------------------------------------- # 

# Mass Action
k = 30
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Simulation

sims = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)

# Plots Generation:

maxval = max([max([max(sims[i][j].I_se[0]+sims[i][j].H_bed[0]+sims[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(alpha[j]))
        #axs[i, j].set_ylim([0,maxval*1.05])
        #axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('MassAction - Total Quarantine')
fig.show()




""" 
 # ---------------------------------------------------------------------- #
 #                                                                        #
 #      Análisis SHFR: Should habe been Hospitalized Fatality Rate        #
 #                                                                        #
 # ---------------------------------------------------------------------- #

    Datos: SHFR
    ejes: VMI totales vs infectados severos 


    Contourplot de SHFR=Muertos acumulados/((Ise+Icr) acumulados) considerando movilidad (alpha) en el eje X, 
     y numero de camas en el eje Y

"""

tstate = ''
initdate = datetime(2020,5,15)


# Parametros del modelo
beta = 0.2 # Tasa de contagio
mu = 0 # Razon E0/I0
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
nm = int(population/100000)
# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0
# UCI/UTI capacity per 100000 persons
Htot_max = 60

step = 2
Htot = list(range(0*nm,Htot_max*nm+step*nm,step*nm))
Htot_per100M = list(range(0,Htot_max+step,step))
# VMI Capacity
Vtot = [i/2 for i in Htot]

# Mobility
step = 0.05
alpha = list(np.arange(0.05,1+step,step))


# ----------------------------------------------------------------- #
#        Plot 1: SHFR Mass Action Dynamics - Total Quarantine       #
# ----------------------------------------------------------------- # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Run  Simulation 
SHFR = np.zeros((len(Htot),len(alpha)))   
sims = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    #aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)                 
 
for i in range(len(Htot)):     
    for j in range(len(alpha)): 
        SHFR[i][j] = sims[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,Htot_per100M,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR Mass Action Dynamics - Total Quarantine')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 


# ------------------------------------------------------------------------ #
#        Plot 2: SHFR Mass Action Dynamics - Dynamic Quarantine (2W)       #
# ------------------------------------------------------------------------ # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
SHFR = np.zeros((len(Htot),len(alpha)))   
sims = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    #aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)                 
 
for i in range(len(Htot)):     
    for j in range(len(alpha)): 
        SHFR[i][j] = sims[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,Htot_per100M,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR Mass Action Dynamics - Dynamic Quarantine (2W)  ')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 



# ------------------------------------------------------------------------ #
#        Plot 3: SHFR Saturated Kynetics Dynamics - Total Quarantine       #
# ------------------------------------------------------------------------ # 
 
 # Mass Action
k = 30
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Run  Simulation 
SHFR = np.zeros((len(Htot),len(alpha)))   
sims = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    #aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)                 
 
for i in range(len(Htot)):     
    for j in range(len(alpha)): 
        SHFR[i][j] = sims[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,Htot_per100M,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR Saturated Kynetics Dynamics - Total Quarantine')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 


# ------------------------------------------------------------------------ #
#        Plot 4: SHFR Saturated Kynetics - Dynamic Quarantine (2W)       #
# ------------------------------------------------------------------------ # 
 
 # Mass Action
k = 30
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
SHFR = np.zeros((len(Htot),len(alpha)))   
sims = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    #aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims.append(aux)                 
 
for i in range(len(Htot)):     
    for j in range(len(alpha)): 
        SHFR[i][j] = sims[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,Htot_per100M,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR Saturated Kynetics - Dynamic Quarantine (2W)')
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
#Grafico se calcula para un valor T de tiempo. Debe ser contourplor con ejes x= alpha, y=kappa,
# y color dado por funcion distancia euclideana entre Infectados kappa =0 y el distinto valor de kappa. 
#Genera contourplots para 4 valores distintos de T: Entre inicio y peak, al peak, entre peak y el fin, y al fin de la dinamica.

"""

tstate = ''
initdate = datetime(2020,5,15)


# Parametros del modelo
beta = 0.2 # Tasa de contagio
mu = 0 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos

tsim = 2000 # Tiempo de simulacion

# Simular
# Activos Iniciales
I_act0 = 10
I_as_prop = 1
I_mi_prop = 0
I_se_prop = 0
I_cr_prop = 0

# Muertos iniciales
dead0 = 0
population = 100000
nm = int(population/100000)
# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0

# UCI/UTI capacity
Htot = 50*nm 
# VMI Capacity
Vtot = Htot/2

# Mobility
step = 0.05
alpha = list(np.arange(0.5,0.75+step,step))

# Saturation Kinetics Factor
step = 2
kmax = 10
k = [0,5,10,20,30,40]#list(np.arange(0,kmax+step,step))


# ------------------------------------------------- #
#        Plot 1: MAD vs SKD -Total Quarantine       #
# ------------------------------------------------- # 

qt = 0
qp = 14
iqt = 0
fqt = tsim
sims = [] 

# Run  Simulation 
sims = [] 
for i in k: 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,Vtot,alpha=j,k=i,qp = qp, qt = qt, iqt = iqt, fqt =fqt) for j in alpha)
    sims.append(aux)


# Euclidian Distance for each alpha for diferent T 

#ED = []
#for i in range(len(k)):
#    aux = []
#    for j in range(len(alpha)):
#        aux.append(LA.norm(sims[i][j].I[0]-sims[0][j].I[0]))
#    ED.append(aux)
        

# Peak size proportion

peak = []
for i in range(len(k)):
    aux = []
    for j in range(len(alpha)):
        aux.append(sims[i][j].peak[0]/sims[0][j].peak[0])
    peak.append(aux)
        

# Contour Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,k,peak) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Peak Size Proportion')
ax.set_xlabel('Mobility')
ax.set_ylabel('Saturation Dynamics Factor')
plt.show() 


# Grid plot
fig, axs = plt.subplots(len(k), len(alpha))
for i in range(len(k)):
    for j in range(len(alpha)):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I[0],label="Infected")  
        axs[i, j].set_title("K: "+str(k[i])+" | Alpha: "+str(alpha[j]))
fig.suptitle('Axes values are scaled individually by default')
#fig.tight_layout()
lines, labels = fig.axes[-1].get_legend_handles_labels()
    
#fig.legend(lines, labels, loc = 'upper center')
fig.legend(lines, labels,loc = 'best')
fig.show()



# SHFR

SHFR = np.zeros((len(k),len(alpha)))
for i in range(len(k)):    
    for j in range(len(alpha)):
        SHFR[i][j] = sims[i][j].SHFR[0] 
        

fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,k,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR')
ax.set_xlabel('Mobility')
ax.set_ylabel('Saturation Dynamics Factor')
plt.show() 


# --------------------------------------------------- #
#        Plot 2: MAD vs SKD -Dynamic Quarantine       #
# --------------------------------------------------- # 


# Run  Simulation 
qt = 1
qp = 14
iqt = 0
fqt = tsim
sims = [] 

for i in k: 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,Vtot,alpha=j,k=i,qp = qp, qt = qt, iqt = iqt, fqt =fqt) for j in alpha)
    sims.append(aux)


# Euclidian Distance for each alpha for diferent T 

ED = []
for i in range(len(k)):
    aux = []
    for j in range(len(alpha)):
        aux.append(LA.norm(sims[i][j].I[0]-sims[0][j].I[0]))
    ED.append(aux)
        

# Peak size proportion

peak = []
for i in range(len(k)):
    aux = []
    for j in range(len(alpha)):
        aux.append(sims[i][j].peak[0]/sims[0][j].peak[0])
    peak.append(aux)
        

# Contour Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,k,peak) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Peak Size Proportion')
ax.set_xlabel('Mobility')
ax.set_ylabel('Saturation Dynamics Factor')
plt.show() 


# Grid plot
fig, axs = plt.subplots(len(k), len(alpha))
for i in range(len(k)):
    for j in range(len(alpha)):
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I[0],label="Infected")  
        axs[i, j].set_title("K: "+str(k[i])+" | Alpha: "+str(alpha[j]))
fig.suptitle('Axes values are scaled individually by default')
#fig.tight_layout()
lines, labels = fig.axes[-1].get_legend_handles_labels()
    
#fig.legend(lines, labels, loc = 'upper center')
fig.legend(lines, labels,loc = 'best')
fig.show()



# SHFR

SHFR = np.zeros((len(k),len(alpha)))
for i in range(len(k)):    
    for j in range(len(alpha)):
        SHFR[i][j] = sims[i][j].SHFR[0] 
        

fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,k,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR')
ax.set_xlabel('Mobility')
ax.set_ylabel('Saturation Dynamics Factor')
plt.show() 




# Comparaciones de todos los graficos anteriores para mass action y saturated kinetics también con distintas cuarentenas
# Fittear datos chilenos 









""" 
Codigo para posible reciclaje



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



fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        #axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_crD_d[0],label="H_cr to D")
        #axs[i, j].plot(sims[i][j].t[0],sims[i][j].VD_d[0],label="V to D")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0]+sims[i][j].I_se[0]+sims[i][j].I_seD_d[0],label="Total Severe")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0],label="Hospitalized")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_seD_d[0],label="I_seD")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        #axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])
fig.suptitle('Axes values are scaled individually by default')
#fig.tight_layout()
lines, labels = fig.axes[-1].get_legend_handles_labels()
    
#fig.legend(lines, labels, loc = 'upper center')
fig.legend(lines, labels,loc = 'best')
fig.show()






sims = []
for i in Vtot:    
    aux = []
    for j in alpha:
        print("alpha")
        # Creación del objeto de simulación 
        simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
        quarantines = [[tsim, 0.85, j, 0.0, 0.0, tsim, 0.0]]
        simulation.inputarray = np.array(quarantines)
        simulation.addquarantine()
        simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,i,R=0,D=0,H_cr = 0)
        simulation.simulate(v=3)
        if i > 0:
            aux.append(simulation)
        else:
            aux.append(simulation)
    sims.append(aux)



num_cores = multiprocessing.cpu_count() 
def ParallelSimulation(Htot,Vtot,alpha): 
    simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop) 
    quarantines = [[tsim, 0.85, alpha, 0.0, 0.0, tsim, 0.0]] 
    simulation.inputarray = np.array(quarantines) 
    simulation.addquarantine() 
    simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
    simulation.simulate(v=3)  
    return simulation #, simulation.SHFR[0]  

"""