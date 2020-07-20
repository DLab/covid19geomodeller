#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# -------------------- #
#                      #
#     SEIRHDV TEsts    #
#                      #
# -------------------- #


Scripts for testing new features and messing around

"""



# -------------------- #
#                      #
#     SEIRHDV Local    #
#                      #
# -------------------- #

from SEIRHDV_local import SEIRHDV_local

















# -------------------- #
#                      #
#     SEIRHDV Remote   #
#                      #
# -------------------- #



from SEIRHDV_remote import SEIRHDV_remote


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
simulation = SEIRHVD_remote(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)

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
