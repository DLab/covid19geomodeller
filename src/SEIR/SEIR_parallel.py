#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
# -------------------- #
#                      #
#    SEIR Parallel     #
#                      #
# -------------------- #

Esta deberia ser la libreria de meta analisis. Hacer una clase por tipo de modelo y agregar los meta-anális necesarios

Estudiar formas más inteligentes de meter argumentos a funciones de python

Agregar iterador sobre k

Luego hay que resolver como generar los meta-plots: Normales, grillas, contour, what-else?  
                                                   


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

from class_SEIR import SEIR



# Get values

class seirMetaAnalysis:
    def __init__(self):
        self.sims = []
        self.num_cores = multiprocessing.cpu_count()
                
    
    def sim_run(self,tsim,alpha,beta,mu,k=0,I=100,I_ac=0,I_d=0,R=0,population=1000000,expinfection = 1, SeroPrevFactor=1,intgr=0): 
        """
            Single SEIR Model simulation
        """
        model = SEIR(tsim,alpha,beta,mu,k=k,I=I,I_ac=I_ac,I_d=I_d,R=R,population=population,expinfection = expinfection, SeroPrevFactor=SeroPrevFactor)      
        if intgr == 0:
            print('Fast Solver')
            model.integr_sci(0,tsim,0.1,False)        
        else:
            print('Robust Solver')            
            model.integr(0,tsim,0.1,False)
        out=model
        return(out)   
    
    def simulate_k(self,tsim,alpha,beta,mu,k=0,I=100,I_ac=0,I_d=0,R=0,population=1000000,expinfection = 1, SeroPrevFactor=1,intgr=0):
        """
            Multi SEIR Model simulation
            
        """        
        self.sims=[]                
        self.sims = (Parallel(n_jobs=self.num_cores, verbose=50)(delayed(self.simulate)(tsim,alpha,beta,mu,k=i,I=I,I_ac=I_ac,I_d=I_d,R=R,population=population,expinfection =expinfection, SeroPrevFactor=SeroPrevFactor,intgr=intgr) for i in k))
        
        self.sims=[]
        if type(alpha) == list:
            self.sims = Parallel(n_jobs=self.num_cores, verbose=50)(delayed(self.sim_run)(tsim,i,beta,mu,k=j,I=I,I_ac=I_ac,I_d=I_d,R=R,population=population,expinfection =expinfection, SeroPrevFactor=SeroPrevFactor,intgr=intgr) for j in k for i in alpha)
        else:
            self.sims = Parallel(n_jobs=self.num_cores, verbose=50)(delayed(self.sim_run)(tsim,alpha,beta,mu,k=j,I=I,I_ac=I_ac,I_d=I_d,R=R,population=population,expinfection =expinfection, SeroPrevFactor=SeroPrevFactor,intgr=intgr) for j in k)
        self.simulated = True
        return(self.sims)



    def simulate_k2(self,tsim,alpha,beta,mu,k=0,I=100,I_ac=0,I_d=0,R=0,population=1000000,expinfection = 1, SeroPrevFactor=1,intgr=0):
        self.sims=[]
        for j in k:
            if type(alpha) == list:
                self.sims.append(Parallel(n_jobs=self.num_cores, verbose=50)(delayed(self.sim_run)(tsim,i,beta,mu,k=j,I=I,I_ac=I_ac,I_d=I_d,R=R,population=population,expinfection =expinfection, SeroPrevFactor=SeroPrevFactor,intgr=intgr) for i in alpha))
            else:
                self.sims.append(self.sim_run(tsim,alpha,beta,mu,k=j,I=I,I_ac=I_ac,I_d=I_d,R=R,population=population,expinfection =expinfection, SeroPrevFactor=SeroPrevFactor,intgr=intgr))
        self.simulated = True
        return(self.sims)



    def simulate(self,tsim,alpha,beta,mu,k=0,I=100,I_ac=0,I_d=0,R=0,population=1000000,expinfection = 1, SeroPrevFactor=1,intgr=0):
        self.sims=[]
        if type(alpha) == list:
            self.sims = Parallel(n_jobs=self.num_cores, verbose=50)(delayed(self.sim_run)(tsim,i,beta,mu,k=k,I=I,I_ac=I_ac,I_d=I_d,R=R,population=population,expinfection =expinfection, SeroPrevFactor=SeroPrevFactor,intgr=intgr) for i in alpha)
        else:
            self.sims = self.sim_run(tsim,alpha,beta,mu,k=k,I=I,I_ac=I_ac,I_d=I_d,R=R,population=population,expinfection =expinfection, SeroPrevFactor=SeroPrevFactor,intgr=intgr)
        self.simulated = True
        return(self.sims)


"""
    def simulate(self,intgr=0):        
        #params=Parallel(n_jobs=num_cores, verbose=50)(delayed(ref_test.refinepso_all)(Ir,tr,swarmsize=200,maxiter=50,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,1],Q_r=[0,1],obj_func='IN')for i in range(int(rep)))
        self.sims=Parallel(n_jobs=self.num_cores, verbose=50)(delayed(self.sim_run)(self.inputarray[i,0],self.inputarray[i,1],self.inputarray[i,2],self.inputarray[i,3],self.inputarray[i,4],self.inputarray[i,5],self.inputarray[i,6]) for i in range(self.inputarray.shape[0]))
        self.simulated = True
        return(self.sims)
"""