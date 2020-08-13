#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -------------------- #
#                      #
#     SEIRHDV Paper    #
#                      #
# -------------------- #
import sys
from pathlib import Path
sys.path.insert(1, '../SEIRHVD/')
sys.path.insert(1, 'SEIRHVD/')
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
 #      Busqueda de Tipping point de movilidad       #
 #                                                   #
 # ------------------------------------------------- #



"""

# Parametros del modelo
beta = 0.2 # Tasa de contagio
mu = 0 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos

# Simular
tsim = 2000 # Tiempo de simulacion
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
Htot = 50*nm

# VMI Capacity
Vtot = 30*nm

#Movilty
# From 0 to 1 in steps of:
step = 0.01
alpha = list(np.arange(0.3,0.35+step,step))

Nalpha = len(alpha)

# Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

quarantines = [[tsim, 0.85, i, qp, iqt, fqt, qt] for i in alpha] 
simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop, k = k) 
simulation.inputarray = np.array(quarantines) 
simulation.addquarantine() 
simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
simulation.simulate(v=3)  

simulation.plotinfectadosactivos()


print('El limite de movilidad para la fijacion de la cuarentena con beta = 0.2 es de alpha=0.33')




""" 
 # ------------------------------------------------- #
 #                                                   #
 #      Busqueda de movilidad basal para 5 camas     #
 #                                                   #
 # ------------------------------------------------- #


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



print('The Basal movility for 5 beds is 0.3')



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
beta = 0.117 # Tasa de contagio
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
NHtot = len(Htot)
# VMI Capacity
Vtot = 10*nm


#Movilty
# From 0 to 1 in steps of:
alpha = [0.65,0.75,0.85]
Nalpha = len(alpha)



# ---------------------------------------------------- #
#        Plot 1:  Cuarentena dinamica 7 dias           #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 7
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
        #axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0]+sims[i][j].I_se[0]+sims[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0],label="Hospitalized")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims[i][j].t[0],sims[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Bed Analysis - Dynamic Quarantine 7 days')
fig.show()


# ---------------------------------------------------- #
#        Plot 2:  Cuarentena dinamica 14 dias        #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim


# Simulation

sims2 = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims2.append(aux)


# Plots Generation:
maxval = max([max([max(sims2[i][j].I_se[0]+sims2[i][j].H_bed[0]+sims2[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = tsim


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        #axs[i, j].plot(sims[i][j].t[0],sims[i][j].H_bed[0]+sims[i][j].I_se[0]+sims[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims2[i][j].t[0],sims2[i][j].H_bed[0],label="Hospitalized")
        axs[i, j].plot(sims2[i][j].t[0],sims2[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims2[i][j].t[0],sims2[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Bed Analysis - Dynamic Quarantine 14 days')
fig.show()

# ---------------------------------------------------- #
#        Plot 3:  Cuarentena dinamica 21 dias          #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 21
qt = 1
iqt = 0
fqt = tsim


# Simulation
sims3 = []
for i in Htot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(i,Vtot,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims3.append(aux)


# Plots Generation:
maxval = max([max([max(sims3[i][j].I_se[0]+sims3[i][j].H_bed[0]+sims3[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NHtot)])
xlim = tsim


fig, axs = plt.subplots(NHtot, Nalpha)
for i in range(NHtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims3[i][j].t[0],sims3[i][j].H_bed[0]+sims[i][j].I_se[0]+sims[i][j].I_seD_d[0],label="Total Severe") # No me convence este grafico
        axs[i, j].plot(sims3[i][j].t[0],sims3[i][j].H_bed[0],label="Hospitalized")
        #axs[i, j].plot(sims3[i][j].t[0],sims3[i][j].I_se[0],label="I_se")
        axs[i, j].plot(sims3[i][j].t[0],sims3[i][j].I_seD_d[0],label="Severe to Death")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))
        axs[i, j].set_ylim([0,maxval*1.05])
        axs[i, j].set_xlim([0,xlim])

lines, labels = fig.axes[-1].get_legend_handles_labels()    
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Bed Analysis - Dynamic Quarantine 21 days')
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
NVtot = len(Vtot)

#Movilty
# From 0 to 1 in steps of:
alpha = [0.65,0.75,0.85]
Nalpha = len(alpha)



# ---------------------------------------------------- #
#           Plot 1:  Cuarentena Total                  #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 0
qt = 0
iqt = 0
fqt = tsim

# Simulation

sims4 = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims4.append(aux)

# Plots Generation:

maxval = max([max([max(sims4[i][j].I_se[0]+sims4[i][j].H_bed[0]+sims4[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims4[i][j].t[0],sims4[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(alpha[j]))
        #axs[i, j].set_ylim([0,maxval*1.05])
        #axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('VMI Analysis - Total Quarantine')
fig.show()



# ---------------------------------------------------- #
#        Plot 2:  Cuarentena dinamica 7 dias         #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 7
qt = 1
iqt = 0
fqt = tsim

# Simulation

sims5 = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims5.append(aux)

# Plots Generation:

maxval = max([max([max(sims5[i][j].I_se[0]+sims5[i][j].H_bed[0]+sims5[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims5[i][j].t[0],sims5[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(alpha[j]))
        #axs[i, j].set_ylim([0,maxval*1.05])
        #axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('VMI Analysis - Dyanmic Quarantine 7 days')
fig.show()


# ---------------------------------------------------- #
#        Plot 3:  Cuarentena dinamica 14 dias         #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Simulation

sims6 = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims6.append(aux)

# Plots Generation:

maxval = max([max([max(sims6[i][j].I_se[0]+sims6[i][j].H_bed[0]+sims6[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims6[i][j].t[0],sims6[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(alpha[j]))
        #axs[i, j].set_ylim([0,maxval*1.05])
        #axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('VMI Analysis - Dyanmic Quarantine 14 days')
fig.show()


# ---------------------------------------------------- #
#        Plot 4:  Cuarentena dinamica 21  dias         #
# ---------------------------------------------------- # 

# Mass Action
k = 0
# Total constant quarantine
qp = 7
qt = 1
iqt = 0
fqt = tsim

# Simulation

sims7 = []
for i in Vtot:
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims7.append(aux)

# Plots Generation:

maxval = max([max([max(sims7[i][j].I_se[0]+sims7[i][j].H_bed[0]+sims7[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NVtot, Nalpha)

for i in range(NVtot):
    for j in range(Nalpha):
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].I_cr[0],label="Critical Infected")
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].VD_d[0],label="Daily VMI to Death")
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].I_crD_d[0],label="Daily Critical Infected to Death")
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].H_crin[0],label="Critical Hospitalized")
        axs[i, j].plot(sims7[i][j].t[0],sims7[i][j].V[0],label="VMI")        
        axs[i, j].set_title("Vtot: "+str(Vtot[i])+" | Alpha: "+str(alpha[j]))
        #axs[i, j].set_ylim([0,maxval*1.05])
        #axs[i, j].set_xlim([0,xlim])        

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('VMI Analysis - Dyanmic Quarantine 21 days')
fig.show()




""" 
 # ---------------------------------------------------------------------- #
 #                                                                        #
 #                           Análisis SHFR(t)                             #
 #                                                                        #
 # ---------------------------------------------------------------------- #

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

NHtot = len(Htot)

#Movilty
# From 0 to 1 in steps of:
alpha = [0.65,0.75,0.85]
Nalpha = len(alpha)



# ------------------------------------------------------------------------ #
#              Plot 1: SHFR(t) Dynamic Quarantine 7 days                 #
# ------------------------------------------------------------------------ # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 7
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
SHFR = np.zeros((len(Htot),len(alpha)))   
sims9 = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    sims9.append(aux)                 

maxval = max([max([max(sims9[i][j].I_se[0]+sims9[i][j].H_bed[0]+sims9[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(len(Htot)):
    for j in range(len(alpha)):
        axs[i, j].plot(sims9[i][j].t[0],sims9[i][j].SHFR_d[0],label="SHFR")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('SHFR(t) - Dynamic Quarantine 7 days')
fig.show()




# ------------------------------------------------------------------------ #
#              Plot 2: SHFR Dynamic Quarantine 14 days                 #
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
sims10 = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    sims10.append(aux)                 

maxval = max([max([max(sims10[i][j].I_se[0]+sims10[i][j].H_bed[0]+sims10[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(len(Htot)):
    for j in range(len(alpha)):
        axs[i, j].plot(sims10[i][j].t[0],sims10[i][j].SHFR_d[0],label="SHFR")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('SHFR - Dynamic Quarantine 7 days')
fig.show()



# ------------------------------------------------------------------------ #
#              Plot 3: SHFR Dynamic Quarantine 21 days                 #
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
sims11 = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    sims11.append(aux)                 

maxval = max([max([max(sims11[i][j].I_se[0]+sims11[i][j].H_bed[0]+sims11[i][j].I_seD_d[0]) for j in range(Nalpha)]) for i in range(NVtot)])
xlim = tsim

fig, axs = plt.subplots(NHtot, Nalpha)

for i in range(len(Htot)):
    for j in range(len(alpha)):
        axs[i, j].plot(sims11[i][j].t[0],sims11[i][j].SHFR_d[0],label="SHFR")
        axs[i, j].set_title("Htot: "+str(Htot[i])+" | Alpha: "+str(alpha[j]))

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('SHFR - Dynamic Quarantine 7 days')
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



# ------------------------------------------------------------------------ #
#              Plot 1: SHFR Dynamic Quarantine 7 days                 #
# ------------------------------------------------------------------------ # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 7
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
  
sims12 = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    #aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims12.append(aux)                 

SHFR = np.zeros((len(Htot),len(alpha))) 
for i in range(len(Htot)):     
    for j in range(len(alpha)): 
        SHFR[i][j] = sims12[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,Htot_per100M,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR Mass Action Dynamics - Dynamic Quarantine 7 days  ')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 




# ------------------------------------------------------------------------ #
#              Plot 2: SHFR Dynamic Quarantine 14 days                 #
# ------------------------------------------------------------------------ # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 
 
sims13 = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    #aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims13.append(aux)                 
 
SHFR = np.zeros((len(Htot),len(alpha)))  
for i in range(len(Htot)):     
    for j in range(len(alpha)): 
        SHFR[i][j] = sims13[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,Htot_per100M,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR Mass Action Dynamics - Dynamic Quarantine 14 days)  ')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 



# ------------------------------------------------------------------------ #
#              Plot 3: SHFR Dynamic Quarantine 21 days                 #
# ------------------------------------------------------------------------ # 
 
 # Mass Action
k = 0
# Total constant quarantine
qp = 14
qt = 1
iqt = 0
fqt = tsim

# Run  Simulation 

sims14 = [] 
for i in range(len(Htot)): 
    aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot[i],Vtot[i],alpha[j],k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in range(len(alpha)))
    #aux = Parallel(n_jobs=num_cores, verbose=50)(delayed(ParallelSimulation)(Htot,i,j,k=k,qp = qp, qt = qt, iqt = iqt, fqt = fqt) for j in alpha)     
    sims14.append(aux)                 
 

SHFR = np.zeros((len(Htot),len(alpha)))   
for i in range(len(Htot)):     
    for j in range(len(alpha)): 
        SHFR[i][j] = sims14[i][j].SHFR[0] 



# Plot
fig,ax=plt.subplots(1,1)
cp = ax.contourf(alpha,Htot_per100M,SHFR) 
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('SHFR Mass Action Dynamics - Dynamic Quarantine 21 days)  ')
ax.set_xlabel('Mobility')
ax.set_ylabel('Beds per 100.000')
plt.show() 
