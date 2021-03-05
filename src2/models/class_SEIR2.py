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
    def __init__(self,tsim,alpha,beta,mu,sigma = 0.2, gamma = 0.1, k=0,I0=100,I_ac0=0,I_d0=0,R0=0,population=1000000,expinfection = 1, SeroPrevFactor=1,chi = 0,psi = 0,k_I=0,k_R=0,RealIC=None, initdate = None,SimIC=None,I_det_prop=1,testaccuracy=1):        
        self.tsim = tsim
        self.alpha = alpha
        self.beta = beta
        self.mu=1.4
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
        
        if type(chi) == int:
            self.chi = np.poly1d(0)
        else:
            self.chi = chi
        
        if type(psi) == int:
            self.psi = np.poly1d(0)
        else:
            self.psi = psi

        self.I_det_prop = I_det_prop # Detected infected propotion
        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #        
        # dVariable/dt = sum(prob_i/in_time_i*in_State_i,i in in states) - sum(prob_i/out_time_i*out_State_i,i in in states) 
        
        # Susceptibles
        # dS/dt:
        self.dS=lambda t,S,E,I,R: self.chi(t) -self.alpha(t)*self.beta*S*(self.expinfection*(E)+I)/(self.N+self.k_I*I+self.k_R*R)+self.eta*R
        # Exposed
        # dE/dt
        self.dE=lambda t,S,E,I,R: self.alpha(t)*self.beta*S*(self.expinfection*(E)+I)/(self.N+self.k_I*I+self.k_R*R)\
            -self.sigma*E

        # Infected
        # dI_as/dt
        self.dI=lambda t,E,I: self.sigma*E - self.gamma*I - self.testaccuracy*self.psi(t)*I/self.population

        # Recovered
        # dR/dt
        self.dR=lambda t,I,R: self.gamma*I+-self.eta*R

        # Acummulated Infected
        self.dI_ac=lambda t,E: self.sigma*E

        # Daily Infected
        self.dI_d = lambda t,E,I_d: self.sigma*E - I_d 
    
        # Performed exams ac
        self.de =  lambda t: self.psi(t)

        # Detected and removed Infected
        self.de_I = lambda t,I: self.testaccuracy*self.psi(t)*I/self.population


    # ------------------- #
    #  Valores Iniciales  #
    # ------------------- #
            
    def setinitvalues(self):       
        self.I = 0        
        self.muS=self.mu        
        self.R=0.0
        self.mu=1.4
        self.t=400.0
 
        self.SeroPrevFactor = 1
        self.population = 8125072
        
        self.expinfection = 1        
        # Accumulated Infected
        self.I_ac = 0

        # Daily Infected
        self.I_d = 0

        # Saturated Kinetics
        self.k = 0

        self.setrelationalvalues()

    def setrelationalvalues(self):        
           
        # Expuestos
        self.E=self.mu*self.I               
        
        # Valores globales
        self.SeroPrevPop =  self.SeroPrevFactor*self.population
        self.S=self.SeroPrevPop-self.E-self.I - self.R
        self.N=(self.S+self.E+self.I+self.R)        

        #constructor of SEIR class elements, it's initialized when a parameter
        #miminization is performed to adjust the best setting of the actual infected


    def integr(self,t0,T,h,E0init=False):
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
            
            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)

        else:
            return()

        
        def model_SEIR_graph(t,y,ydot):
            
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3])
            ydot[1]=self.dE(t,y[0],y[1],y[2],y[3])
            ydot[2]=self.dI(t,y[1],y[2])
            ydot[3]=self.dR(t,y[2],y[3])
            
            ydot[4]=self.dI_ac(t,y[1])
            ydot[5]=self.dI_d(t,y[1],y[5])

            ydot[6]=self.de(t)
            ydot[7]=self.de_I(t,y[2])

            
        initcond = np.array([S0,E0,I0,R0,I_ac0,I_d0,e0,e_I0])                                


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

        #Cálculo de la fecha del Peak  
        self.peakindex = np.where(self.I==max(self.I))[0][0]
        self.peak = max(self.I)
        self.peak_t = self.t[self.peakindex]
        if self.initdate:
            self.peak_date = self.initdate+timedelta(days=round(self.peak_t)) 

        # Prevalence: 
        self.prevalence_total = self.I_ac/self.population
        self.prevalence_susc = [self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]
        self.prevalence_det = [self.I_det*self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]                         
               
        return(sol)

    def integr_sci(self,t0,T,h,E0init=False):
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

            self.t=np.arange(t0,T+h,h)

        else:
            return()
            

        
        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3])
            ydot[1]=self.dE(t,y[0],y[1],y[2],y[3])
            ydot[2]=self.dI(t,y[1],y[2])
            ydot[3]=self.dR(t,y[2],y[3])
            
            ydot[4]=self.dI_ac(t,y[1])
            ydot[5]=self.dI_d(t,y[1],y[5])

            ydot[6]=self.de(t)
            ydot[7]=self.de_I(t,y[2])

            return(ydot)
        initcond = np.array([S0,E0,I0,R0,I_ac0,I_d0,e0,e_I0])  

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








