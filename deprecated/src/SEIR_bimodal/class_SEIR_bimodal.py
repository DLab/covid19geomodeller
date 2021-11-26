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

"""
To do:
  - Create reports function inside simSAEIRHVD class


SEIRHVD Implementation
Instructions: 
    Init a simSEIRHVD objecting giving the simulation condictions:
        - tsim: Simulation time
        - max_mov:
        - rem_mov:
        - qp:
        - iqt:
        - fqt:
        - movfunct: 

"""


class SEIRbimodal:  
    def __init__(self,tsim,alpha_s,alpha_r,beta,mu,sigma = 0.2,gamma = 0.1,p_r=0.5,k=0,I0=100,I_ac0=0,I_d0=0,R=0,population=1000000,expinfection = 1, SeroPrevFactor=1):
        
        self.tsim = tsim
        self.alpha_s = alpha_s
        self.alpha_r = alpha_r
        self.beta = beta
        self.mu=mu
        self.k = k
        self.p_r = p_r                

        self.sigma = sigma # Incubation rate (period^⁻1)
        self.gamma = gamma # Recovery rate (period^⁻1)        

        self.Is = I0*(1-p_r)
        self.Ir = I0*p_r
        self.Is_ac = I_ac0*(1-p_r)
        self.Ir_ac = I_ac0*p_r
        self.Is_d = I_d0*(1-p_r)
        self.Ir_d = I_d0*p_r

        self.R = R 
 
        self.SeroPrevFactor = SeroPrevFactor
        self.population = population        
        self.expinfection = expinfection         
        
        self.t=0

        # Valores iniciales
        # Expuestos
        self.Es = self.mu*self.Is
        self.Er = self.mu*self.Ir
        
        # Valores globales
        self.N =  self.SeroPrevFactor*self.population        
        self.Ss = self.N*(1-self.p_r) - self.Es-self.Is - self.R*(1-self.p_r)
        self.Sr = self.N*(self.p_r) - self.Er-self.Ir - self.R*(self.p_r)

                     
        
        
        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #        
        # dVariable/dt = sum(prob_i/in_time_i*in_State_i,i in in states) - sum(prob_i/out_time_i*out_State_i,i in in states) 
        
        # Susceptibles
        # dSs/dt:
        self.dSs=lambda t,Ss,Es,Er,Is,Ir: -self.alpha_s(t)*self.beta*Ss*(self.expinfection*(Es+Er)+Is+Ir)/(self.N+self.k*(Is+Ir))
        # dSr/dt:
        self.dSr=lambda t,Sr,Es,Er,Is,Ir: -self.alpha_s(t)*self.beta*Sr*(self.expinfection*(Es+Er)+Is+Ir)/(self.N+self.k*(Is+Ir)) - self.beta*self.alpha_r(t)*Sr*Ir/(self.p_r*self.N)
        # Exposed
        # dEs/dt
        self.dEs=lambda t,Ss,Es,Er,Is,Ir: self.alpha_s(t)*self.beta*Ss*(self.expinfection*(Es+Er)+Is+Ir)/(self.N+self.k*(Is+Ir)) - self.sigma*Es
        # dEr/dt
        self.dEr=lambda t,Sr,Es,Er,Is,Ir: self.alpha_s(t)*self.beta*Sr*(self.expinfection*(Es+Er)+Is+Ir)/(self.N+self.k*(Is+Ir)) + self.beta*self.alpha_r(t)*Sr*Ir/(self.p_r*self.N) \
            - self.sigma*Er 
            

        # Infected        
        self.dIs=lambda t,Es,Is: self.sigma*Es - self.gamma*Is
        self.dIr=lambda t,Er,Ir: self.sigma*Er - self.gamma*Ir

        # Recovered
        # dR/dt
        self.dR=lambda t,Is,Ir: self.gamma*(Is+Ir)

        # Acumulated Infected:
        self.dIs_ac=lambda t,Es: self.sigma*Es
        self.dIr_ac=lambda t,Er: self.sigma*Er

        # Daily Infected:
        self.dIs_d=lambda t,Es,Is_d: self.sigma*Es - Is_d
        self.dIr_d=lambda t,Er,Ir_d: self.sigma*Er - Ir_d
    
    def globalvars(self):        
        self.S = self.Ss+self.Sr
        self.E = self.Er + self.Es
        self.I = self.Is+self.Ir
        
    def integr(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.
        print('Import scikits-odes')
        from scikits.odes.odeint import odeint


        if(not isinstance(self.Ss, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                Es0=self.mu*(self.Is)
                Er0=self.mu*(self.Ir)
            else:
                Es0=self.Es
                Er0=self.Er
                
            Ss0 = self.Ss
            Sr0 = self.Sr

            Is0 = self.Is
            Ir0 = self.Ir

            R0 = self.R            
            
            Is_ac0 = self.Is_ac
            Ir_ac0 = self.Ir_ac
            
            Is_d0 = self.Is_d
            Ir_d0 = self.Ir_d

            self.t=np.arange(t0,T+h,h)
            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition
            Ss0 = self.Ss[idx]
            Sr0 = self.Sr[idx]

            Es0=self.Es[idx]
            Er0=self.Er[idx]

            Is0 = self.Is[idx]
            Ir0 = self.Ir[idx]

            R0 = self.R[idx]
            
            Is_ac0 = self.Is_ac[idx]
            Ir_ac0 = self.Ir_ac[idx]
            
            Is_d0 = self.Is_d[idx]
            Ir_d0 = self.Ir_d[idx]
                        
            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)

        else:
            return()

        
        def model_SEIR_graph(t,y,ydot):
            #Ss
            ydot[0]=self.dSs(t,y[0],y[2],y[3],y[4],y[5])
            #Sr
            ydot[1]=self.dSr(t,y[1],y[2],y[3],y[4],y[5])
            #Es
            ydot[2]=self.dEs(t,y[0],y[2],y[3],y[4],y[5])
            #Er
            ydot[3]=self.dEr(t,y[1],y[2],y[3],y[4],y[5])
            #Is
            ydot[4]=self.dIs(t,y[2],y[4])
            #Ir
            ydot[5]=self.dIr(t,y[3],y[5])
            #R
            ydot[6]=self.dR(t,y[4],y[5])
            #Is_ac
            ydot[7]=self.dIs_ac(t,y[2])
            #Ir_ac
            ydot[8]=self.dIr_ac(t,y[3])
            #Is_d
            ydot[9]=self.dIs_d(t,y[2],y[9])
            #Ir_d
            ydot[10]=self.dIr_d(t,y[3],y[10])            

        initcond = np.array([Ss0,Sr0,Es0,Er0,Is0,Ir0,R0,Is_ac0,Ir_ac0,Is_d0,Ir_d0])



        sol = odeint(model_SEIR_graph, self.t, initcond,method='admo')
        
        self.t=sol.values.t         
        self.Ss=sol.values.y[:,0]
        self.Sr=sol.values.y[:,1]
        self.Es=sol.values.y[:,2]
        self.Er=sol.values.y[:,3]
        self.Is=sol.values.y[:,4]
        self.Ir=sol.values.y[:,5]
        self.R=sol.values.y[:,6]
        self.Is_ac=sol.values.y[:,7]
        self.Ir_ac=sol.values.y[:,8]
        self.Is_d=sol.values.y[:,9]
        self.Ir_d=sol.values.y[:,10]

        self.globalvars()              
        return(sol)

    def integr_sci(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.

        if(not isinstance(self.Ss, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                Es0=self.mu*(self.Is)
                Er0=self.mu*(self.Ir)
            else:
                Es0=self.Es
                Er0=self.Er
                
            Ss0 = self.Ss
            Sr0 = self.Sr

            Is0 = self.Is
            Ir0 = self.Ir

            R0 = self.R            
            
            Is_ac0 = self.Is_ac
            Ir_ac0 = self.Ir_ac
            
            Is_d0 = self.Is_d
            Ir_d0 = self.Ir_d

            self.t=np.arange(t0,T+h,h)

            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition
            Ss0 = self.Ss[idx]
            Sr0 = self.Sr[idx]

            Es0=self.Es[idx] # Maybe I need to take idx out
            Er0=self.Er[idx] # Maybe I need to take idx out

            Is0 = self.Is[idx]
            Ir0 = self.Ir[idx]

            R0 = self.R[idx]
            
            Is_ac0 = self.Is_ac[idx]
            Ir_ac0 = self.Ir_ac[idx]
            
            Is_d0 = self.Is_d[idx]
            Ir_d0 = self.Ir_d[idx]
                        
            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)


        else:
            return()
            

        
        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            #Ss
            ydot[0]=self.dSs(t,y[0],y[2],y[3],y[4],y[5])
            #Sr
            ydot[1]=self.dSr(t,y[1],y[2],y[3],y[4],y[5])
            #Es
            ydot[2]=self.dEs(t,y[0],y[2],y[3],y[4],y[5])
            #Er
            ydot[3]=self.dEr(t,y[1],y[2],y[3],y[4],y[5])
            #Is
            ydot[4]=self.dIs(t,y[2],y[4])
            #Ir
            ydot[5]=self.dIr(t,y[3],y[5])
            #R
            ydot[6]=self.dR(t,y[4],y[5])
            #Is_ac
            ydot[7]=self.dIs_ac(t,y[2])
            #Ir_ac
            ydot[8]=self.dIr_ac(t,y[3])
            #Is_d
            ydot[9]=self.dIs_d(t,y[2],y[9])
            #Ir_d
            ydot[10]=self.dIr_d(t,y[3],y[10]) 
                                          
            return(ydot)

        initcond = np.array([Ss0,Sr0,Es0,Er0,Is0,Ir0,R0,Is_ac0,Ir_ac0,Is_d0,Ir_d0])  

        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA')
        
        self.t=sol.t 
        
        self.Ss=sol.y[0,:]
        self.Sr=sol.y[1,:]
        self.Es=sol.y[2,:]
        self.Er=sol.y[3,:]
        self.Is=sol.y[4,:]
        self.Ir=sol.y[5,:]
        self.R=sol.y[6,:]
        self.Is_ac=sol.y[7,:]
        self.Ir_ac=sol.y[8,:]
        self.Is_d=sol.y[9,:]
        self.Ir_d=sol.y[10,:]

        self.globalvars()
        return(sol)


