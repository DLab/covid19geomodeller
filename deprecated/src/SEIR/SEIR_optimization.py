#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linalg as LA
import pygmo as pg

import pandas as pd
from time import time

from class_SEIR import SEIR
from Quarantine import Quarantine

"""
    Testing pygmo applied into SEIR Class
"""
# To do: 
#  - Revisar una mejor construccion del objeto SEIRModel, usaria el objeto seir que hicimos al principio
class SEIRModel_pygmo:
    def __init__(self,Ir,tr,S0,I0,R0,h,mov,qp,movefunct,bounds):
        self.Ir = Ir
        self.tr = tr
        self.S0 = S0
        self.I0 = I0
        self.R0 = R0
        self.h = h
        self.mov = mov
        self.qp = qp
        self.movefunct = movefunct
        self.bounds = bounds
    def fitness(self,x):        
        self.E0=x[3]*self.I0
        sol=pd.DataFrame(SDSEIR.intger(self.S0,self.E0,self.I0,self.R0,min(self.tr),max(self.tr),self.h,x[0],x[1],x[2],self.mov,self.qp,self.tr[-1],self.movefunct))
        idx=np.searchsorted(sol.t,self.tr)
        res = LA.norm(self.Ir-sol.I[idx])        
        return([res])

    def get_bounds(self):
        return(self.bounds)

    def set_bounds(self,bounds):
        self.bounds = bounds
        return(self.bounds)

    def get_name(self):
        return "SEIR Unsectorial"



# PSO Normal
algo = pg.algorithm(pg.pso(gen = 20))
pop = pg.population(prob,50)
t0 = time()
pop = algo.evolve(pop)
t1 = time()
print('Optimization takes %f seconds' %(t1-t0))
print(pop.champion_f)
print(pop.champion_x)





class SEIRModel:
    def __init__(self,Ir,tr,I_ac_r,I_ac_tr,tsim,alpha,population,bounds):
        self.Ir = Ir
        self.tr = tr
        self.I_ac_r = I_ac_r
        self.I_ac_tr = I_ac_tr
        self.tsim = tsim        
        self.alpha = alpha
        self.population = population        
        self.bounds = bounds
    def fitness(self,x):        
        sol = SEIR(tsim=self.tsim,alpha=self.alpha,beta=x[0],mu=x[1],k=x[2],I=self.Ir[0],I_ac=self.I_ac_r[0],population=self.population)
        sol.integr_sci(0,tsim,0.1)                
        idx=np.searchsorted(sol.t,self.tr)
        res = LA.norm(self.Ir-sol.I[idx])
        return([res])

    def get_bounds(self):
        return(self.bounds)

    def set_bounds(self,bounds):
        self.bounds = bounds
        return(self.bounds)




from datetime import datetime
from importdata import ImportData as importadata
# Import data: 
tstate = '13'
initdate = datetime(2020,5,15)

Ir,tr,Ir_dates = importdata.importActiveInfected(tstate = tstate, initdate = initdate)
I_ac_r,I_ac_tr,I_ac_dates = importdata.importAcumulatedInfected(tstate = tstate, initdate = initdate)
population = ImportData.importPopulation(tstate = tstate)



tsim = 1000
mob = 0.6
iqt = 0
alpha = Quarantine(mob,iqt=iqt).alpha

# Params to find
# beta,mu,k
lb=[0.01,0.1, 0]
ub=[   1,  2,30]
bounds = [lb,ub]

opti = SEIRModel(Ir=Ir,tr=tr,I_ac_r=I_ac_r,I_ac_tr=I_ac_tr,tsim=tsim,alpha=alpha,population=population,bounds=bounds)

algo = pg.algorithm(pg.pso(gen = 20))
pop = pg.population(opti,50)
t0 = time()
pop = algo.evolve(pop)
t1 = time()
print('Optimization takes %f seconds' %(t1-t0))
print(pop.champion_f)
print(pop.champion_x)