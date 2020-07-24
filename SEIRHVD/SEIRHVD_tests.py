#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# -------------------- #
#                      #
#     SEIRHDV TEsts    #
#                      #
# -------------------- #


Scripts for testing new features and messing around.
This commands must be run from SEIRHVD directory

"""



# -------------------- #
#                      #
#     SEIRHDV Local    #
#                      #
# -------------------- #

from SEIRHVD_local import SEIRHVD_local
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

# ------------------------------------------- #
#        Región Metropolitana      #
# ------------------------------------------- #


# Región Por CUT 
tstate = '13'
# Fecha Inicial 
initdate = datetime(2020,4,13)

# Fecha de inicio de cuarentena
qit =  datetime(2020,5,15)
qit = (qit - initdate).days

qfd = datetime(2020,7,24)
qft = (qfd - initdate).days

# Parametros del modelo
beta = 0.24#25#15#15#2#117
mu = 0.1#5 
ScaleFactor = 3
SeroPrevFactor = 1#0.22 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos
k=35 # Factor de Saturación Cinética

# Tiempo de simulacion
tsim = 1000 

# Creación del objeto de simulación 
simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,k=k)

# Creación del vector de cuarentenas
# [Tsim, max_mov,rem_mov,quarantine period, quarantine initial time, quarantine final time, quarantine type]
quarantines = [[500.0, 0.85, 0.6, 0.0, qit,qft, 0.0],
               [500.0, 0.85, 0.65, 0.0,qit,qft, 0.0],
               [500.0, 0.85, 0.7, 0.0, qit,qft, 0.0]]
               #[500.0, 0.85, 0.75, 0.0,qit,qft, 0.0],
               #[500.0, 0.85, 0.4, 0.0,qit, qft, 0.0]]
simulation.inputarray = np.array(quarantines) # This will change during next update
simulation.addquarantine()

simulation.simulate(v=3)

# Plots de ejemplo:
simulation.plotinfectadosactivos(scalefactor=True,ylim=200000) 11

simulation.plotinfectadosactivos() 
simulation.plotventiladores()
simulation.plotfallecidosacumulados()


# ------------------..........--------------- #
#        Simulación Datos Artificiales        #
# ------------------------------------------- #

tstate = ''
# Fecha Inicial
initdate = datetime(2020,5,15)
# Parametros del modelo
beta = 0.2 # Tasa de contagio
mu = 0.6 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos
tsim = 500 # Tiempo de simulacion

# Creación del objeto de simulación 
simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)


quarantines = [[500.0, 0.85, 0.7, 0.0, 0.0, 500.0, 0.0]]
simulation.inputarray = np.array(quarantines)
simulation.addquarantine()

# Valores iniciales
I_act0 = 100 
dead = 0
population = 10000
H0 = 1
V0 = 1
# Capacidades hospitalarias
Htot = 1000
Vtot = 1000

simulation.initialvalues(I_act0,dead,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)

# Simular
simulation.simulate()

# Plots de ejemplo:
simulation.plotfallecidosdiarios()
simulation.plotventiladores()


# QA:

simulation.I_as_ac_prop[0][-1] +simulation.I_mi_ac_prop[0][-1] +simulation.I_se_ac_prop[0][-1] +simulation.I_cr_ac_prop[0][-1] +simulation.I_act0/simulation.Iac[0][-1] 






# -------------------- #
#                      #
#     SEIRHDV Remote   #
#                      #
# -------------------- #
# Not ready yet


from SEIRHVD_remote import SEIRHVD_remote


# ------------------..........--------------- #
#        Ingreso de Parámetros Generales      #
# ------------------------------------------- #
# Región Por CUT 
tstate = '13'

# Fecha Inicial
initdate = datetime(2020,5,15)

# Parametros del modelo
beta = 0.117 # Tasa de contagio
mu = 0.6 # Razon E0/I0
ScaleFactor = 1.9 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 0.5 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos

tsim = 500 # Tiempo de simulacion

# Creación del objeto de simulación 
simulation = SEIRHVD_DA(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)

# Creación de escenarios
# Opcion 1: Escenarios por defecto 
simulation.defaultescenarios()

# Opción 2: Escenarios personalizados
# [tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct]

inputarray = [[500.0, 0.85, 0.6, 0.0, 0.0, 500.0, 0.0],
              [500.0, 0.85, 0.65, 0.0, 0.0, 500.0, 0.0],
              [500.0, 0.85, 0.7, 0.0, 0.0, 500.0, 0.0]]
simulation.inputarray = np.array(inputarray)
simulation.addscenario()


# Simular
simulation.simulate()

# Análisis de Resultados
# funcion genérica: 
simulation.pĺotvariable(enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False)

simulation.plotactivos()
# Infectados Activos  
simulation.tabladedatos(inicio = datetime(2020,5,15), fin = datetime(2020,6,30),variables =['I_cum','I_act','D','L'], path=''))
