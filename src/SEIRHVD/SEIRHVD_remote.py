#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import dill as pickle
import json
import requests
from datetime import datetime

import SEIRHVD_importdata
import SEIRHVD_tables
import SEIRHVD_plots
import SEIRHVD_quarantine
import SEIRHVD_vars
from datetime import datetime
from datetime import timedelta

"""
# ------------------------------------------- #
#                                             #
#            SEIRHDV Remote Simulation        #
#                                             #
# ------------------------------------------- #
    
Remote simmulation using ZEUS Cluster   
Connection to Dlab's VPN is necessary in order for to work remotely

To Do:
- Send data as dill/pickle

"""

class SEIRHVD_remote(SEIRHVD_tables.SEIRHVD_tables,SEIRHVD_plots.SEIRHVD_plots,SEIRHVD_importdata.SEIRHVD_importdata,SEIRHVD_vars.SEIRHVD_vars,SEIRHVD_quarantine.SEIRHVD_quarantine):        
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

        if tstate:            
            self.inputdata = True
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
           

    # ----------------------- #
    #     Run simmulation     #
    # ----------------------- #
    
    def simulate(self,v=3,intgr=0):
        if not self.inputdata:
            return('Set initnial values before simulating')
        endpoint = 'http://192.168.2.248:5003/SEIRHVDsimulate'
        auxinputarray = [list(self.inputarray[i]) for i in range(self.numescenarios)]
        data = {
        'state': str(self.tstate),
        'beta': str(self.beta),
        'mu': str(self.mu),        
        'tsim': str(self.tsim),
        'initdate': self.initdate.strftime('%Y/%m/%d'),
        'ScaleFactor': str(self.ScaleFactor),
        'SeroPrevFactor': str(self.SeroPrevFactor),
        'inputarray': str(auxinputarray),
        'version':v,
        'intgr':intgr,
        'B':self.B,
        'D':self.D,
        'population':self.population,
        'I_act0':self.I_act0,
        'H0':self.H0,
        'V':self.V,
        'Htot':self.Htot,
        'Vtot':self.Vtot,
        'H_cr':self.H_cr,
        'R':self.R}
       
        pickle.dumps(data)
        r = requests.post(url = endpoint, data = data)
        self.sims = pickle.loads(r.content)
        self.importdata()                 
        self.localvar()      
        return