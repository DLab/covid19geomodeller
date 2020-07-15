#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ----------------- #
#                   #
#    SEIR Model     #
#                   #
# ----------------- #



from SEIRHVD_local import SEIRHVD_local
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing



tstate = ''
initdate = datetime(2020,5,15)

# ------------------- #
#        Plot 1       #
# ------------------- # 

# Camas H totales vs infectados severos vs HcrtoD




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

population = 1000000
# Initial Hospitalized
H0 = 0
# Initial VMI 
V0 = 0
# UCI/UTI capacity per 1000 persons
nm = int(population/100000)
#step = 6


# Hospital capacity
Htot = 30*nm

# VMI Capacity
Vtot = 10*nm



# Quarantines
alpha = 0.6
quarantines = [[tsim, 0.85, alpha, 0.0, 0.0, tsim, 0.0]]



# Simulation

# Kinetic Saturation
k=0

simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate, k = k)
simulation.inputarray = np.array(quarantines)
simulation.addquarantine()
simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
simulation.simulate(v=0)

# Kinetic Saturation
k=0.1


simulation2 = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate, k = k)
simulation2.inputarray = np.array(quarantines)
simulation2.addquarantine()
simulation2.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
simulation2.simulate(v=0)