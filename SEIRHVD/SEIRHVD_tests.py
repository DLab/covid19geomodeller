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


# Parametros del modelo
beta = 0.117
mu = 0.6 
ScaleFactor = 1.9
SeroPrevFactor = 0.22 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos


# Tiempo de simulacion
tsim = 500 

# Creación del objeto de simulación 
simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)

# Creación del vector de cuarentenas
# [Tsim, max_mov,rem_mov,quarantine period, quarantine initial time, quarantine final time, quarantine type]
quarantines = [[500.0, 0.85, 0.6, 0.0, qit, 500.0, 0.0],
              [500.0, 0.85, 0.65, 0.0,qit, 500.0, 0.0],
              [500.0, 0.85, 0.7, 0.0, qit, 500.0, 0.0],
              [500.0, 0.85, 0.75, 0.0,qit, 500.0, 0.0],
              [500.0, 0.85, 0.4, 0.0,qit, 500.0, 0.0]]
simulation.inputarray = np.array(quarantines) # This will change during next update
simulation.addquarantine()

simulation.simulate(v=3)

# Plots de ejemplo:
simulation.plotinfectadosacumulados() 
simulation.plotventiladores()
simulation.plotfallecidosacumulados()


# ------------------..........--------------- #
#        Simulación Datos Artificiales        #
# ------------------------------------------- #

tstate = ''
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
simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)


quarantines = [[500.0, 0.85, 0.65, 0.0, 0.0, 500.0, 0.0]]
simulation.inputarray = np.array(quarantines)
simulation.addquarantine()

# Valores iniciales
I_act0 = 1000 
dead = 10 
population = 100000
H0 = 10
V0 = 1
# Capacidades hospitalarias
Htot = 20
Vtot = 10

simulation.initialvalues(I_act0,dead,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)

# Simular
simulation.simulate()

# Plots de ejemplo:
simulation.plotfallecidosdiarios()
simulation.plotventiladores()







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
