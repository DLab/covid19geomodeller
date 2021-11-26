import numpy as np
from numpy import linalg as LA 
import matplotlib.pyplot as plt

from pathlib import Path
sys.path.insert(1, '../src/SEIR/')
sys.path.insert(1, '../SEIR/')
sys.path.insert(1, '../src/utils/')
sys.path.insert(1, '../utils/')

from class_SEIR_bimodal import SEIRbimodal
from class_SEIR import SEIR
from Quarantine import Quarantine


"""
SEIR bimodal QA 

SEIR bimodal has 2 types of population interacting in a single geographical unit

"""

# Parameters

beta = 0.3
mu = 1 
# Incubation rate
gamma = 0.2
# Recovery rate
sigma = 0.1
# Risk population proportion
p_r = 0.4
# Total population
population = 1000000

# Saturation Dynamics:
k = 0

# Mobility
mob_s = 0.8
mob_r = 0.8

alpha_s = Quarantine(mob_s).alpha
alpha_r = Quarantine(mob_r).alpha

# simulation time
tsim = 1000

# Initial conditions
I0 = 200
I_ac0 = 0
I_d0 = 0

# Do exposed infect?
expinfection = 0

# Create model and integrate
seirbimodal = SEIRbimodal(tsim,alpha_s,alpha_r,beta,mu,sigma = sigma,gamma = gamma,p_r=p_r,k=k,I0=I0,I_ac0=I_ac0,I_d0=I_d0,population=population,expinfection = expinfection )

seirbimodal.integr_sci(0,tsim,0.1)


# Create local vars:

t = seirbimodal.t
Ss = seirbimodal.Ss
Sr = seirbimodal.Sr
S = Ss+Sr
Es = model.Es
Er = model.Er
E = Er + Es
Is = model.Is
Ir = model.Ir
I = Is+Ir
R = model.R

# Plot Total People
plt.plot(seirbimodal.t,seirbimodal.S+seirbimodal.E+seirbimodal.I+seirbimodal.R,label='Simulation total people')
plt.hlines(population,0,tsim,label='Total People') 
plt.title('Population QA')
plt.legend(loc=0)
plt.show()

# Plot SEIR
plt.plot(t,S,label='S')
plt.plot(t,E,label='E')
plt.plot(t,I,label='I')
plt.plot(t,R,label='R')
plt.title('SEIR')
plt.legend(loc=0)
plt.show()

# Plot Infected
plt.plot(t,Is,label='Is')
plt.plot(t,Ir,label='Ir')
plt.title('Active Infected')
plt.legend(loc=0)
plt.show()
