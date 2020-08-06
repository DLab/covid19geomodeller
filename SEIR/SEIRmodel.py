#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import class_SEIR as SEIR
import SEIR_plots
import SEIR_vars
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

"""

class SEIRmodel(SEIR_plots.SEIR_plots,SEIR_vars.SEIR_vars,SEIRHVD_quarantine.SEIRHVD_quarantine):        
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
        #if tstate:
        #    self.importdata()           
        #    self.I_act0 = self.ScaleFactor*self.Ir[0]            
        #    self.inputdata = True
        #    Hcmodel = np.poly1d(np.polyfit(self.sochimi_tr, self.Hr_tot, bedmodelorder)) 
        #    Vcmodel = np.poly1d(np.polyfit(self.sochimi_tr, self.Vr_tot, bedmodelorder))
        #    Vmax = np.mean(self.Vr_tot[-7:])*1.01# 1500 Vcmodel(tsat)
        #    Hmax = np.mean(self.Hr_tot[-7:])*1.01# self.Hcmodel(tsat)            
        #    try:
        #        vents = [Vcmodel(t) for t in range(self.tsim)]          
        #        tsat = int(np.where(np.array(vents)>=Vmax)[0][0])
        #    except:
        #        tsat = self.sochimi_tr[-1]+7            
        #    self.Htot=lambda t: Hcmodel(t)*(1-expit(t-tsat)) + expit(t-tsat)*Hmax  #1997.0        
        #    self.Vtot=lambda t: Vcmodel(t)*(1-expit(t-tsat)) + expit(t-tsat)*Vmax    
        #    self.H0 = self.Hr[0] # Hospitalizados totales dia 0  
        #    self.V = self.Vr[0]   # Ventilados al dia de inicio
        #    self.H_cr = 0 #Hospitalizados a la espera de un ventilador día 0
        #    self.B = self.Br[0]  # Muertos acumulados al dia de inicio
        #    self.D = self.Br[1]-self.Br[0]  # Muertos en el dia de inicio
        #    self.R  = 0
        #    
        #    self.realdata = True
        if True:
            self.inputdata = False
            print('Set initial values')
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
        


    # -------------- #
    #    Simulate    #
    # -------------- #

    def simulate(self,intgr=0):
        if not self.inputdata:
            return('Set the initial values before running the simulation!')

        print('SEIR Model')
        model = SEIR.simSEIRHVD(beta = self.beta, mu = self.mu, inputarray= self.inputarray, B=self.B,D=self.D,V=self.V,I_act0=self.I_act0,R=self.R,Htot=self.Htot,Vtot=self.Vtot,H_cr=self.H_cr,H0=self.H0,expinfection=self.expinfection, SeroPrevFactor= self.SeroPrevFactor, population = self.population,intgr=intgr,I_as_prop = self.I_as_prop, I_mi_prop = self.I_mi_prop,I_se_prop = self.I_se_prop,I_cr_prop = self.I_cr_prop,k=self.k)
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