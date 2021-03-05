#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIRHVD Model
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
from joblib import Parallel, delayed
from scipy import signal
import pandas as pd
from numpy import linalg as LA 
import multiprocessing  
#import SEIR_plots

from datetime import datetime
from datetime import timedelta

"""
SEIRQ Model Implementation

"""

class SEIR:  
    def __init__(self,tsim,beta,gamma = 0.1, I0=100,I_ac0=0,I_d0=0,R0=0,population=1000000, initdate = None):        
        self.tsim = tsim
        
        self.beta = beta

        self.gamma = gamma # Recovery rate
        
        self.I = I0
        self.I_ac = I_ac0
        self.I_d = I_d0

        self.R = R0
 
        self.population = population        
        
        self.t=0
       
        # Valores globales
        self.N = self.population
        self.S = self.N- self.I - self.R   

        self.initdate = None 

        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #        
        # dVariable/dt = sum(prob_i/in_time_i*in_State_i,i in in states) - sum(prob_i/out_time_i*out_State_i,i in in states) 
        
        # Susceptibles
        # dS/dt:
        self.dS=lambda t,S,I: -self.beta*S*I/self.N
    
        # Infected
        # dI_as/dt
        self.dI=lambda t,S,I: -self.beta*S*I/self.N - self.gamma*I 

        # Recovered
        # dR/dt
        self.dR=lambda t,I: self.gamma*I

        # Acummulated Infected
        self.dI_ac=lambda t,S,I: self.beta*S*I/self.N

        # Daily Infected
        self.dI_d = lambda t,S,I,I_d: self.beta*S*I/self.N - I_d 
    

    def integr_sci(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.
       
        if(not isinstance(self.S, np.ndarray)):

            S0=self.S
            I0=self.I
            R0=self.R            
            
            I_ac0=self.I_ac
            I_d0=self.I_d

            self.t=np.arange(t0,T+h,h)
            
        #elif((min(self.t)<=t0) & (t0<=max(self.t))):
        #    #Condition over exiting time in already initialized object

        #    #Search fot initial time
        #    idx=np.searchsorted(self.t,t0)

        #    #set initial condition

        #    S0=self.S[idx]
        #    E0=self.E
        #    
        #    I0=self.I[idx]
        #    R0=self.R[idx]                        
        #    I_ac0=self.I_ac[idx]
        #    I_d0=self.I_d[idx]                       

        #    e0 = self.e[idx]
        #    e_I0 = self.e_I[idx]

        #    self.t=np.arange(t0,T+h,h)

        else:
            return()
            

        
        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dS(t,y[0],y[1])
            ydot[1]=self.dI(t,y[0],y[1])
            ydot[2]=self.dR(t,y[1])            
            ydot[3]=self.dI_ac(t,y[0],y[1])
            ydot[4]=self.dI_d(t,y[0],y[1],y[4])

            return(ydot)
        initcond = np.array([S0,I0,R0,I_ac0,I_d0])  

        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA')
        
        self.t=sol.t 

         
        self.S=sol.y[0,:]
        self.I=sol.y[1,:]
        self.R=sol.y[2,:]
        self.I_ac=sol.y[3,:]
        self.I_d=sol.y[4,:]

        #Cálculo de la fecha del Peak  
        self.peakindex = np.where(self.I==max(self.I))[0][0]
        self.peak = max(self.I)
        self.peak_t = self.t[self.peakindex]

        # Prevalence: 
        self.prevalence_total = self.I_ac/self.population
        return(sol)








