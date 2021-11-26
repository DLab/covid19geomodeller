#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from SEIRHVD_local import SEIRHVD_local
from datetime import datetime
import numpy as np

# ----------------------------- #
#        Model Parameters       #
# ----------------------------- #

# Region code (Chile)
tstate = '13'
# Simulation Initial date
initdate = datetime(2020,4,13)

# Epidemiological Paramters
beta = 0.2#117
mu = 0.9 # Initial Exposed/Infected Proportion 
ScaleFactor = 1.9 # Under-reported cases proportion
SeroPrevFactor = 1 # Sero Prevalence Correction Factor
expinfection = 1 # Exposed contagion efectivity 
k=10 # Kinetic Saturation Factor 

# Simulation time
tsim = 500 

# Simulation Object Construction
simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,k=k)



# Quarantines:
  
# Quaratine Initial date
qid =  datetime(2020,5,15)
qit = (qid - initdate).days
# Quaratine Final date
qfd = datetime(2020,7,24)
qft = (qfd - initdate).days
# [Tsim, max_mov,rem_mov,quarantine period, quarantine initial time, quarantine final time, quarantine type]
quarantines = [[500.0, 0.85, 0.6, 0.0, qit,qft, 0.0],
               [500.0, 0.85, 0.65, 0.0,qit,qft, 0.0],
               [500.0, 0.85, 0.7, 0.0, qit,qft, 0.0]]

simulation.inputarray = np.array(quarantines) # This will change during next update
simulation.addquarantine()


# Simulate
simulation.simulate(v=3)

# Results
simulation.plotinfectadosactivos(days = -1, scalefactor=True)
simulation.plotinfectadosacumulados(days = -1, scalefactor=True) 
simulation.plotventiladores(days = -1)
simulation.plotfallecidosacumulados(days = -1)
simulation.plotcuarentenas()


# Dynamical Quarantines
# Creación del objeto de simulación 
simulation2 = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate,k=k)

# Creación del vector de cuarentenas
# [Tsim, max_mov,rem_mov,quarantine period, quarantine initial time, quarantine final time, quarantine type]
qit = (qid - initdate).days
qft = (qfd - initdate).days

quarantines = [[500.0, 0.85, 0.65, 0.0, qit,qft, 0],
               [500.0, 0.85, 0.65, 7,qit,qft, 1],
               [500.0, 0.85, 0.65, 14, qit,qft, 1]]

simulation2.inputarray = np.array(quarantines) # This will change during next update
simulation2.addquarantine()

# Simulate
simulation2.simulate(v=3)

# Results
simulation2.plotinfectadosactivos(days = -1, scalefactor=True)
simulation2.plotinfectadosacumulados(days = -1, scalefactor=True) 
simulation2.plotventiladores(days = -1)
simulation2.plotfallecidosacumulados(days = -1)
simulation2.plotcuarentenas()
