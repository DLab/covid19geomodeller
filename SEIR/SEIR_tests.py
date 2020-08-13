#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ----------------- #
#                   #
#    SEIR Model     #
#                   #
# ----------------- #

import sys
from pathlib import Path
sys.path.insert(1, '../src/')
sys.path.insert(1, 'src/')
from SEIRmodel import SEIRmodel
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
#%matplotlib tk
from joblib import Parallel, delayed
import multiprocessing

#from SEIRHVD_local import SEIRHVD_local



tstate = ''
initdate = datetime(2020,5,15)

# ------------------- #
#        Plot 1       #
# ------------------- # 

beta = 0.12 # Contagion rate
mu = 0.5 # E0/I0 initial rate

k=40 # Kinetic Saturation: 0 for mass action mixing

ScaleFactor = 1 # Scale Factor: Number of real infected over reported 
SeroPrevFactor = 1 # Sero Prevalence Factor: Adjust the proportion of the population that enters the virus dynamics
expinfection = 1 # Exposed contagion rate compared to the infected (0 the don't infect, 1 the infect in the same rate as the infected )

# Simulation time
tsim = 1000
# Population in case we don't pick a specific place
if not tstate:
    population = 1000000
# Initial Active Infected 
I_act0 = 100


# Quarantines
max_mob = 0.8 # Maximum mobility
# Total quarantine
s1 = [tsim,max_mob,0.65,0,0,tsim,0]
s2 = [tsim,max_mob,0.5,0,0,tsim,0]
# Dynamic quarantine
s3 = [tsim,max_mob,0.3,14,0,tsim,1]
s4 = [tsim,max_mob,0.5,14,0,tsim,1]

# Define one quarantine array for each different quarantine remanent mobility
quarantines = [s1,s3]





simulation = SEIRmodel(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,tsim = tsim,tstate=tstate, k = k,initdate=initdate)
simulation.inputarray = np.array(quarantines)
simulation.addquarantine()
if not tstate:
    simulation.initialvalues(I_act0,population,R=0)


simulation.simulate()











# Kinetic Saturation
k=0

simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate, k = k,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
simulation.inputarray = np.array(quarantines)
simulation.addquarantine()
simulation.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
simulation.simulate(v=0)

simulation.plotseird()
simulation.peak
simulation.peak_t

k=0.1

simulation2 = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate, k = k,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
simulation2.inputarray = np.array(quarantines)
simulation2.addquarantine()
simulation2.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
simulation2.simulate(v=0)

simulation2.plotseird()
simulation2.peak
simulation2.peak_t

k=1

simulation3 = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate, k = k,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
simulation3.inputarray = np.array(quarantines)
simulation3.addquarantine()
simulation3.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
simulation3.simulate(v=0)

simulation3.plotseird()
simulation3.peak
simulation3.peak_t

k=10

simulation4 = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate, k = k,I_as_prop = I_as_prop, I_mi_prop = I_mi_prop,I_se_prop = I_se_prop,I_cr_prop = I_cr_prop)
simulation4.inputarray = np.array(quarantines)
simulation4.addquarantine()
simulation4.initialvalues(I_act0,dead0,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
simulation4.simulate(v=0)

simulation4.plotseird()
simulation4.peak
simulation4.peak_t

# Plots:

plt.plot(simulation.t[0],simulation.I[0],label="k=0")
plt.plot(simulation2.t[0],simulation2.I[0],label="k=0.1")
plt.plot(simulation3.t[0],simulation3.I[0],label="k=1")
plt.plot(simulation4.t[0],simulation4.I[0],label="k=10")
plt.legend(loc=0)
plt.ylabel('Infected')
plt.xlabel('Days')
plt.xlim(0,500)
plt.show()











"""

 SEIR 2 tests

"""



import sys
from pathlib import Path
sys.path.insert(1, '../SEIR/')
sys.path.insert(1, 'SEIR/')

import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing

from class_SEIR2 import SEIR
from Quarantine import Quarantine


tsim = 1000
max_mob = 0.85
rem_mob = 0.4 
qp = 0
iqt = 20
alpha = Quarantine.alphafunct(rem_mob,iqt = iqt)

beta = 0.117
mu = 1.5
k = 40
I_ac = 100
I = 100
population = 1000000

model = SEIR(tsim,alpha,beta,mu,k=0,I=I,I_ac=I_ac,population=population)