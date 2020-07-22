#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import class_SEIRHUVD3 as SD3
import class_SEIRHUVD2 as SD2
import class_SEIR as SEIR
import SEIRHVD_importdata
import SEIRHVD_tables
import SEIRHVD_plots
import SEIRHVD_quarantine
import SEIRHVD_vars
from datetime import datetime
from datetime import timedelta
import numpy as np
from scipy.special import expit

"""
# ------------------------------------------- #
#                                             #
#            SEIRHDV Local Simulation         #
#                                             #
# ------------------------------------------- #

Requires scikitis.odes library

"""

class SEIRHVD_local(SEIRHVD_tables.SEIRHVD_tables,SEIRHVD_plots.SEIRHVD_plots,SEIRHVD_importdata.SEIRHVD_importdata,SEIRHVD_vars.SEIRHVD_vars,SEIRHVD_quarantine.SEIRHVD_quarantine):
        
    def __init__(self,beta,mu,ScaleFactor=1,SeroPrevFactor=1,expinfection=1,initdate = datetime(2020,5,15), tsim = 500,tstate='',bedmodelorder=2, k = 0,I_as_prop = 0.35, I_mi_prop = 0.63,I_cr_prop = 0.007,I_se_prop = 0.013):
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
        self.I_as_prop = I_as_prop
        self.I_mi_prop = I_mi_prop
        self.I_cr_prop = I_cr_prop
        self.I_se_prop = I_se_prop          
        
        # Import real data
        if tstate:
            self.importdata()           
            self.I_act0 = self.ScaleFactor*self.Ir[0]            
            self.inputdata = True
            Hcmodel = np.poly1d(np.polyfit(self.sochimi_tr, self.Hr_tot, bedmodelorder)) 
            Vcmodel = np.poly1d(np.polyfit(self.sochimi_tr, self.Vr_tot, bedmodelorder))
            Vmax = np.mean(self.Vr_tot[-7:])*1.01# 1500 Vcmodel(tsat)
            Hmax = np.mean(self.Hr_tot[-7:])*1.01# self.Hcmodel(tsat)            
            try:
                vents = [Vcmodel(t) for t in range(self.tsim)]          
                tsat = int(np.where(np.array(vents)>=Vmax)[0][0])
            except:
                tsat = self.sochimi_tr[-1]+7            
            self.Htot=lambda t: Hcmodel(t)*(1-expit(t-tsat)) + expit(t-tsat)*Hmax  #1997.0        
            self.Vtot=lambda t: Vcmodel(t)*(1-expit(t-tsat)) + expit(t-tsat)*Vmax    
            self.H0 = self.Hr[0] # Hospitalizados totales dia 0  
            self.V = self.Vr[0]   # Ventilados al dia de inicio
            self.H_cr = 0 #Hospitalizados a la espera de un ventilador día 0
            self.B = self.Br[0]  # Muertos acumulados al dia de inicio
            self.D = self.Br[1]-self.Br[0]  # Muertos en el dia de inicio
            self.R  = 0
            
            self.realdata = True
        else:
            self.inputdata = False
            print('Set initial values')
            self.realdata = False
        return

    def initialvalues(self,I_act0,dead,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0):
        self.B=dead
        self.D = D
        self.population = population
        self.I_act0 = I_act0
        self.H0=H0
        self.V=V0        
        self.Htot = np.poly1d(Htot) 
        self.Vtot = np.poly1d(Vtot)
        self.H_cr = H_cr
        self.R  = R
        self.inputdata = True
        


    # -------------- #
    #    Simulate    #
    # -------------- #

    def simulate(self,v=3,intgr=0):
        if not self.inputdata:
            return('Set the initial values before running the simulation!')

        if v==2:
            model = SD2.simSEIRHVD(beta = self.beta, mu = self.mu, inputarray= self.inputarray, B=self.B,D=self.D,V=self.V,I_act0=self.I_act0,R=self.R,Htot=self.Htot,Vtot=self.Vtot,H_cr=self.H_cr,H0=self.H0,expinfection=self.expinfection, SeroPrevFactor= self.SeroPrevFactor, population = self.population,intgr=intgr)
        elif v==3:
            model = SD3.simSEIRHVD(beta = self.beta, mu = self.mu, inputarray= self.inputarray, B=self.B,D=self.D,V=self.V,I_act0=self.I_act0,R=self.R,Htot=self.Htot,Vtot=self.Vtot,H_cr=self.H_cr,H0=self.H0,expinfection=self.expinfection, SeroPrevFactor= self.SeroPrevFactor, population = self.population,intgr=intgr,I_as_prop = self.I_as_prop, I_mi_prop = self.I_mi_prop,I_se_prop = self.I_se_prop,I_cr_prop = self.I_cr_prop,k=self.k)
        elif v==0:
            print('SEIR Model')
            model = SEIR.simSEIRHVD(beta = self.beta, mu = self.mu, inputarray= self.inputarray, B=self.B,D=self.D,V=self.V,I_act0=self.I_act0,R=self.R,Htot=self.Htot,Vtot=self.Vtot,H_cr=self.H_cr,H0=self.H0,expinfection=self.expinfection, SeroPrevFactor= self.SeroPrevFactor, population = self.population,intgr=intgr,I_as_prop = self.I_as_prop, I_mi_prop = self.I_mi_prop,I_se_prop = self.I_se_prop,I_cr_prop = self.I_cr_prop,k=self.k)
        else:
            raise('Version Error')
        self.sims = model.simulate()
        self.localvar()
        return