from SEIRHVD_local import SEIRHVD_local
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
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
simulation.plotinfectadosactivos()