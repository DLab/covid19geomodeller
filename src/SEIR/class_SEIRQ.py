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
    def __init__(self,tsim,beta,alpha=1,mu=1.4,sigma = 0.2, gamma = 0.1, k=0,I0=100,I_ac0=0,I_d0=0,R0=0,population=1000000,expinfection = 0, SeroPrevFactor=1,chi = 0,psi = 0,k_I=0,k_R=0,RealIC=None, initdate = None,SimIC=None,I_det_prop=1,testaccuracy=1,lambda_Q=1,T_T = 1,T_Q = 14,lambda_Tr=0):        
        self.tsim = tsim        
        self.beta = beta
        self.mu=mu
        self.k_I = k_I
        self.k_R = k_R

        self.sigma = sigma # Incubation rate
        self.gamma = gamma # Recovery rate
        self.eta = 0.0 # Immunity loss rate

        self.I = I0
        self.I_ac = I_ac0
        self.I_d = I_d0

        self.R = R0
 
        self.SeroPrevFactor = SeroPrevFactor
        self.population = population        
        self.expinfection = expinfection         
        
        self.t=0

        # Expuestos
        self.E = self.mu*self.I               
        
        # Valores globales
        self.N =  self.SeroPrevFactor*self.population
        self.S = self.N- self.E-self.I - self.R   

        self.numescenarios = 1
        self.initdate = None 

        self.testaccuracy = testaccuracy
        
        if type(alpha) == int:
            self.alpha = np.poly1d(alpha)
        else:
            self.alpha = alpha


        if type(chi) == int:
            self.chi = np.poly1d(0)
        else:
            self.chi = chi
        

        # Exams and individual quarantines
        if type(psi) == int:
            self.psi = np.poly1d(0)
        else:
            self.psi = psi

        self.I_det_prop = I_det_prop # Detected infected propotion
        self.lambda_Q = lambda_Q # Efectively Quarantined proportion (1: every detected is quarantined, 0: None of them) 
        self.T_T = T_T # Time between test and Results -> Quarantine
        self.T_Q = T_Q # Quarantine Time
        self.lambda_Tr = lambda_Tr # Tracing effect. This shows the additional inffected located due to tracing.

        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #        
        # dVariable/dt = sum(prob_i/in_time_i*in_State_i,i in in states) - sum(prob_i/out_time_i*out_State_i,i in in states) 
        
        # Susceptibles
        # dS/dt:
        self.dS=lambda t,S,E,I,R,I_T: self.chi(t) -self.alpha(t)*self.beta*S*(self.expinfection*(E)+I+I_T)/(self.N+self.k_I*(I+I_T)+self.k_R*R)+self.eta*R
        # Exposed
        # dE/dt
        self.dE=lambda t,S,E,I,R,I_T: self.alpha(t)*self.beta*S*(self.expinfection*(E)+I+I_T)/(self.N+self.k_I*(I+I_T)+self.k_R*R)\
            -self.sigma*E

        # Infected
        # dI_as/dt
        self.dI=lambda t,E,I: self.sigma*E - self.gamma*I - self.testaccuracy*self.lambda_Q*self.psi(t)*I/self.population*(1+self.lambda_Tr)

        # Recovered
        # dR/dt
        self.dR=lambda t,I,R,Q: self.gamma*I+-self.eta*R + Q/self.T_Q


        # Acummulated Infected
        self.dI_ac=lambda t,E: self.sigma*E

        # Daily Infected
        self.dI_d = lambda t,E,I_d: self.sigma*E - I_d 
    
        # Performed exams ac
        self.de =  lambda t: self.psi(t)

        # Detected and removed Infected
        self.de_I = lambda t,I: self.testaccuracy*self.psi(t)*I/self.population


        # Infected that have been tested but not quarantined yet 
        self.dI_T=lambda t,I,I_T: self.testaccuracy*self.lambda_Q*self.psi(t)*I/self.population*(1+self.lambda_Tr) - I_T/self.T_T

        # Quarantined
        self.dQ= lambda t,I_T,Q: I_T/self.T_T - Q/self.T_Q


    def integr(self,t0=0,T=500,h=0.01,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.
        print('Import scikits-odes')
        from scikits.odes.odeint import odeint


        if(not isinstance(self.S, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                E0=self.mu*(self.I)                
            else:
                E0=self.E
                
            S0=self.S
            I0=self.I
            R0=self.R            
            
            I_ac0=self.I_ac
            I_d0=self.I_d

            e0 = 0
            e_I0 = 0

            I_T0 = 0 
            Q0 = 0 

            self.t=np.arange(t0,T+h,h)
            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            S0=self.S[idx]
            E0=self.E
            
            I0=self.I[idx]
            R0=self.R[idx]                        
            I_ac0=self.I_ac[idx]
            I_d0=self.I_d[idx]

            e0 = self.e[idx]
            e_I0 = self.e_I[idx]

            I_T0 = self.I_T[idx]
            Q0 = self.Q[idx]
            
            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)

        else:
            return()

        
        def model_SEIR_graph(t,y,ydot):
            
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3],y[8])
            ydot[1]=self.dE(t,y[0],y[1],y[2],y[3],y[8])
            ydot[2]=self.dI(t,y[1],y[2])
            ydot[3]=self.dR(t,y[2],y[3],y[9])
            
            ydot[4]=self.dI_ac(t,y[1])
            ydot[5]=self.dI_d(t,y[1],y[5])

            ydot[6]=self.de(t)
            ydot[7]=self.de_I(t,y[2])

            ydot[8]=self.dI_T(t,y[2],y[8])
            ydot[9]=self.dQ(t,y[8],y[9])


            
        initcond = np.array([S0,E0,I0,R0,I_ac0,I_d0,e0,e_I0,I_T0,Q0])


        sol = odeint(model_SEIR_graph, self.t, initcond,method='admo')
        
        self.t=sol.values.t 
        
        self.S=sol.values.y[:,0]
        self.E=sol.values.y[:,1]        
        self.I=sol.values.y[:,2]
        self.R=sol.values.y[:,3]
        self.I_ac=sol.values.y[:,4]
        self.I_d=sol.values.y[:,5]
        self.e=sol.values.y[:,6]
        self.e_I=sol.values.y[:,7]

        self.I_T=sol.values.y[:,8]
        self.Q=sol.values.y[:,9]

        #Cálculo de la fecha del Peak  
        self.peakindex = np.where(self.I==max(self.I))[0][0]
        self.peak = max(self.I)
        self.peak_t = self.t[self.peakindex]
        if self.initdate:
            self.peak_date = self.initdate+timedelta(days=round(self.peak_t)) 

        # Detected cases
        self.I_det = self.I_det_prop*self.I
        self.I_ac_det = self.I_det_prop*self.I_ac
        self.I__d_det = self.I_det_prop*self.I_d


        # Prevalence: 
        self.prevalence_total = self.I_ac/self.population
        self.prevalence_susc = [self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]
        self.prevalence_det = [self.I_det*self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]                         
               
        return(sol)

    def integr_sci(self,t0=0,T=500,h=0.01,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.
       
        if(not isinstance(self.S, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                E0=self.mu*(self.I)                
            else:
                E0=self.E
                
            S0=self.S
            I0=self.I
            R0=self.R            
            
            I_ac0=self.I_ac
            I_d0=self.I_d

            e0 = 0
            e_I0 = 0
            I_T0 = 0 
            Q0 = 0            

            self.t=np.arange(t0,T+h,h)
            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            S0=self.S[idx]
            E0=self.E
            
            I0=self.I[idx]
            R0=self.R[idx]                        
            I_ac0=self.I_ac[idx]
            I_d0=self.I_d[idx]                       

            e0 = self.e[idx]
            e_I0 = self.e_I[idx]

            I_T0 = self.I_T[idx]
            Q0 = self.Q[idx]            

            self.t=np.arange(t0,T+h,h)

        else:
            return()
            

        
        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3],y[8])
            ydot[1]=self.dE(t,y[0],y[1],y[2],y[3],y[8])
            ydot[2]=self.dI(t,y[1],y[2])
            ydot[3]=self.dR(t,y[2],y[3],y[9])
            
            ydot[4]=self.dI_ac(t,y[1])
            ydot[5]=self.dI_d(t,y[1],y[5])

            ydot[6]=self.de(t)
            ydot[7]=self.de_I(t,y[2])

            ydot[8]=self.dI_T(t,y[2],y[8])
            ydot[9]=self.dQ(t,y[8],y[9])            

            return(ydot)
        initcond = np.array([S0,E0,I0,R0,I_ac0,I_d0,e0,e_I0,I_T0,Q0])  

        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA')
        
        self.t=sol.t 

         
        self.S=sol.y[0,:]
        self.E=sol.y[1,:]
        self.I=sol.y[2,:]
        self.R=sol.y[3,:]
        self.I_ac=sol.y[4,:]
        self.I_d=sol.y[5,:]

        self.e=sol.y[6,:]
        self.e_I=sol.y[7,:]

        self.I_T=sol.y[8,:]
        self.Q=sol.y[9,:]

        #Cálculo de la fecha del Peak  
        self.peakindex = np.where(self.I==max(self.I))[0][0]
        self.peak = max(self.I)
        self.peak_t = self.t[self.peakindex]
        if self.initdate:
            self.peak_date = self.initdate+timedelta(days=round(self.peak_t))    

        # Detected cases
        self.I_det = self.I_det_prop*self.I
        self.I_ac_det = self.I_det_prop*self.I_ac
        self.I__d_det = self.I_det_prop*self.I_d

        # Prevalence: 
        self.prevalence_total = self.I_ac/self.population
        self.prevalence_susc = [self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]
        self.prevalence_det = [self.I_det_prop*self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]
        return(sol)








