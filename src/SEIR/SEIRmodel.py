#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import class_SEIR as SEIR
import SEIR_plots
import SEIR_vars
import SEIR_importdata
from datetime import datetime
from datetime import timedelta
import numpy as np
from scipy.special import expit
import SEIRHVD_quarantine

"""
# ------------------------------------------- #
#                                             #
#              SEIR Simulation                #
#                                             #
# ------------------------------------------- #

Old version that uses parallelization inside class_SEIR.py

"""

class SEIRmodel(SEIR_plots.SEIR_plots,SEIR_vars.SEIR_vars,SEIRHVD_quarantine.SEIRHVD_quarantine,SEIR_importdata.SEIR_importdata):        
    def __init__(self,beta,mu,ScaleFactor=1,SeroPrevFactor=1,expinfection=1,tsim = 500,tstate='',bedmodelorder=2, k = 0,initdate = datetime(2020,5,15)):
        self.beta = beta
        self.mu = mu
        self.ScaleFactor = ScaleFactor
        self.SeroPrevFactor = SeroPrevFactor
        self.expinfection = expinfection
        self.tstate = tstate
        self.initdate = initdate
        self.tsim = tsim
        self.May15 = (datetime(2020,5,15)-initdate).days
        self.k = k
        self.I_as_prop = 1
        self.I_mi_prop = 0
        self.I_cr_prop = 0
        self.I_se_prop = 0          
        
        # Import real data
        if tstate:
            self.importdata()           
            self.I_act0 = self.ScaleFactor*self.Ir[0]            
            self.inputdata = True
            self.R  = 0            
            self.realdata = True
        else:
            self.inputdata = False
            print('You must set initial values')
            self.realdata = False
        return

    def initialvalues(self,I_act0,population,R=0):
        self.B=0
        self.D = 0
        self.population = population
        self.I_act0 = I_act0
        self.H0=0
        self.V=0        
        self.Htot = np.poly1d(1000) 
        self.Vtot = np.poly1d(1000)
        self.H_cr = 0
        self.R  = R
        self.inputdata = True
        print('Initial values setted')
        


    # -------------- #
    #    Simulate    #
    # -------------- #

    def simulate(self,intgr=0):
        if not self.inputdata:
            return('Set the initial values before running the simulation!')

        print('SEIR Model')
        self.Htot = np.poly1d(1000) 
        self.Vtot = np.poly1d(1000)        
        model = SEIR.simSEIRHVD(beta = self.beta, mu = self.mu, inputarray= self.inputarray, I_act0=self.I_act0,R=self.R,expinfection=self.expinfection, SeroPrevFactor= self.SeroPrevFactor, population = self.population,intgr=intgr,k=self.k,Htot=self.Htot,Vtot=self.Vtot)
        self.sims = model.simulate()
        self.localvar()
        return

    """
    def ParallelSimulation(alpha,k=0, qp = 0, qt = 0,iqt = 0,fqt = 500): 
        simulation = SEIRmodel(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,tsim = tsim,tstate='', k = k)
        quarantines = [[tsim, 0.85, alpha, qp, iqt, fqt, qt]] 
        simulation.inputarray = np.array(quarantines) 
        simulation.addquarantine() 
        simulation.initialvalues(I_act0,population,R=0)
        simulation.simulate()  
        return simulation

    """