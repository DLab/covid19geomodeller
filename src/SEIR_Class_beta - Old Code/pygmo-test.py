#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linalg as LA
import pygmo as pg
import Single_dist_ref_SEIR as SDSEIR
import pandas as pd
from time import time


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

# PSO Mejorada
algo = pg.algorithm(pg.pso_gen(gen = 50,memory = True, variant = 6))
pop = pg.population(prob,50)
t0 = time()
pop = algo.evolve(pop)
t1 = time()
print('Optimization takes %f seconds' %(t1-t0))
print(pop.champion_f)
print(pop.champion_x)


# Estudio
algo = pg.algorithm(pg.de(gen = 50))
pop = pg.population(prob,50)
t0 = time()
pop = algo.evolve(pop)
t1 = time()
print('Optimization takes %f seconds' %(t1-t0))
print(pop.champion_f)
print(pop.champion_x)

# Comparacion metodos: 
# PSO se demora menos que el bee_colony y llega a resultados simulares