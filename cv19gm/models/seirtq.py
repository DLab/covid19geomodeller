#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIRTQ Model
"""

import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
from datetime import timedelta
import utils.cv19files as cv19files


class SEIRTQ:
    """
        SEIRTQ model object:
        Construction:
            SEIRTQ(self, config = None, inputdata=None)

    """
    def __init__(self, config = None, inputdata=None,verbose = False, **kwargs):
        self.compartmentalmodel = "SEIRTQ"
        self.verbose = verbose
        
        if not config:
            print('Missing configuration file, using default')
    
        # Load Parameters
        if self.verbose:
            print('Loading configuration file')          
        cv19files.loadconfig(self,config,inputdata,**kwargs)
        
        if self.verbose:
            print('Initializing parameters and variables')
        self.set_initial_values()
        if self.verbose:
            print('Building equations')          
        self.set_equations()

        self.solved = False
        if self.verbose:
            print('SEIR object created')
  
    def set_initial_values(self):
        # Underrreporting 
        # Active infected
        if hasattr(self,'I_det'):
            self.I = self.I_det/self.pI_det
        else:
            self.I_det = self.I*self.pI_det

        # New daily Infected
        if hasattr(self,'I_d_det'):
            self.I_d = self.I_d_det/self.pI_det
        else:
            self.I_d_det = self.I_d*self.pI_det


        # Accumulated Infected
        if hasattr(self,'I_ac_det'):
            self.I_ac = self.I_ac_det/self.pI_det
        else:
            self.I_ac_det = self.I_ac*self.pI_det

        
        # Exposed
        #if not self.Einit:
        if not hasattr(self,'E'):
            self.E = self.mu*self.I
        elif not self.E:
            self.E = self.mu*self.I

        if not hasattr(self,'E_d'):
            self.E_d=self.mu*self.I_d
        elif not self.E_d:
            self.E_d=self.mu*self.I_d

        if not hasattr(self,'E_ac'):    
            self.E_ac=self.mu*self.I_ac
        elif not self.E_ac:
            self.E_ac=self.mu*self.I_ac
       
        # Valores globales
        if not hasattr(self,'popfraction'):
            self.popfraction = 1
        # Unitialized states
        if not hasattr(self,'R_d'):
            self.R_d = 0        
        if not hasattr(self,'T'):
            self.T = 0        
        if not hasattr(self,'T_d'):
            self.T_d = 0     
        if not hasattr(self,'Q'):
            self.Q = 0     

        # Population        
        self.N = self.popfraction*self.population
        self.S = self.N - self.E - self.I - self.R - self.T - self.Q

        # External flux functions: 
        if not hasattr(self,'S_f'):
            self.S_f = lambda t:0
        if not hasattr(self,'E_f'):
            self.E_f = lambda t:0
        if not hasattr(self,'I_f'):
            self.I_f = lambda t:0 
        if not hasattr(self,'R_f'):
            self.R_f = lambda t:0
        
    def set_equations(self):
        """        
        Sets Diferential Equations
        """
        # --------------------------- #
        #        Susceptibles         #
        # --------------------------- #
        # 0) dS/dt:
        self.dS=lambda t,S,I,R,T,N: self.S_f(t) - self.alpha(t)*self.beta(t)*S*(I+T)/N + self.rR_S(t)*R
        
        # --------------------------- #
        #           Exposed           #
        # --------------------------- #     
        # 1) dE/dt
        self.dE = lambda t,S,E,I,T,N: self.E_f(t) + self.alpha(t)*self.beta(t)*S*(I+T)/N - E/self.tE_I(t)
 
        # 2) Daily dE/dt
        self.dE_d = lambda t,S,E_d,I,T,N: self.E_f(t) + self.alpha(t)*self.beta(t)*S*(I+T)/N - E_d

        # --------------------------- #
        #           Infected          #
        # --------------------------- #                
        # 3) Active
        self.dI=lambda t,E,I,N: self.I_f(t) + E/self.tE_I(t) - I/self.tI_R(t) - self.k_Ex(t)*self.k_Eacc(t)*self.k_Q(t)*(1+self.k_Tr(t))*I/N
        
        # 4) New Daily
        self.dI_d = lambda t,E,I_d: self.I_f(t) + E/self.tE_I(t) - I_d

        # --------------------------- #
        #         Recovered           #
        # --------------------------- #  
        # 5) Total recovered
        self.dR=lambda t,I,R,Q: self.R_f(t) + I/self.tI_R(t) + Q/self.tQ_R(t) - self.rR_S(t)*R

        # 6) Recovered per day
        self.dR_d=lambda t,I,R_d,Q: self.R_f(t) + I/self.tI_R(t) + Q/self.tQ_R(t) - R_d

        # --------------------------------- #
        #      Testing and Quarantines      #
        # --------------------------------- #  
        # 7) Tested Infected:
        self.dT=lambda t,I,T,N: self.T_f(t) + self.k_Ex(t)*self.k_Eacc(t)*self.k_Q(t)*(1+self.k_Tr(t))*I/N - T/self.tT_Q(t)

        # 8) Daily Tested Infected:
        self.dT_d=lambda t,I,T_d,N: self.T_f(t) + self.k_Ex(t)*self.k_Eacc(t)*self.k_Q(t)*(1+self.k_Tr(t))*I/N - T_d

        # 9) Quarantined:
        self.dQ=lambda t,T,Q:  T/self.tT_Q(t) - Q/self.tQ_R(t)

        # 10) External Flux: // Revisar esto! 
        self.dN = lambda t: self.S_f(t) + self.E_f(t) + self.I_f(t) + self.R_f(t) + self.T_f(t)

    def run(self,t0=0,T=None,h=0.01):
        self.solve(t0=t0,T=T,h=h)

    def solve(self,t0=0,T=None,h=0.01,method='LSODA'):
        """
        Solves ODEs using scipy.integrate
        Args:
            t0 (int, optional): Initial time. Defaults to 0.
            T ([type], optional): Endtime. Defaults to time given when building the object
            h (float, optional): Time step. Defaults to 0.01.
        """

        if T is None:
            T = self.tsim

        # Check if we already simulated the array
        if self.solved:
            #print('Already solved')
            return()
        

        self.t=np.arange(t0,T+h,h)
        initcond = np.array([self.S,self.E,self.E_d,self.I,self.I_d,self.R,self.R_d,self.T,self.T_d,self.Q,self.N]) # [S0,E0,E_d0,I0,I_d0,R0,R_d0,Flux0]
        
        sol = solve_ivp(self.model_equations,(t0,T), initcond,method=method,t_eval=list(range(t0,T)))
        
        self.sol = sol
        self.t=sol.t 
        
        self.S=sol.y[0,:]
        self.E=sol.y[1,:]
        self.E_d=sol.y[2,:]
        self.I=sol.y[3,:]
        self.I_d=sol.y[4,:]
        self.R=sol.y[5,:]
        self.R_d=sol.y[6,:]
        self.T=sol.y[7,:]
        self.T_d=sol.y[8,:]
        self.Q=sol.y[9,:]
        self.N=sol.y[10,:]

        self.E_ac = np.cumsum(self.E_d)
        self.I_ac = np.cumsum(self.I_d) + self.I_ac # second term is the initial condition
        self.R_ac = np.cumsum(self.R_d)      
        self.T_ac = np.cumsum(self.T_d) 

        self.I_det = self.I*self.pI_det
        self.I_d_det = self.I_d*self.pI_det
        self.I_ac_det = self.I_ac*self.pI_det

        self.analytics()
        self.df_build()
        self.solved = True
        return 

    def model_equations(self,t,y):
        ydot=np.zeros(len(y))
        ydot[0]=self.dS(t,y[0],y[3],y[5],y[7],y[10])
        ydot[1]=self.dE(t,y[0],y[1],y[3],y[7],y[10])
        ydot[2]=self.dE_d(t,y[0],y[2],y[3],y[7],y[10])
        ydot[3]=self.dI(t,y[1],y[3],y[10])
        ydot[4]=self.dI_d(t,y[1],y[4])
        ydot[5]=self.dR(t,y[3],y[5],y[9])
        ydot[6]=self.dR_d(t,y[3],y[6],y[9])
        ydot[7]=self.dT(t,y[3],y[7],y[10])
        ydot[8]=self.dT_d(t,y[3],y[8],y[10])
        ydot[9]=self.dQ(t,y[7],y[9])
        ydot[10]=self.dN(t)                         
        return(ydot)

    def analytics(self):
        """
        Perform simulation analytics after running it.
        It calculates peaks, prevalence, and will include R(t). 
        """
        #Cálculo de la fecha del Peak  
        self.peakindex = np.where(self.I==max(self.I))[0][0]
        self.peak = max(self.I)
        self.peak_t = self.t[self.peakindex]
        if self.initdate:
            self.dates = [self.initdate+timedelta(int(self.t[i])) for i in range(len(self.t))]
            self.peak_date = self.initdate+timedelta(days=round(self.peak_t)) 
        else:
            self.dates = [None for i in range(len(self.t))]
            self.peak_date = None            

        # Prevalence: 
        self.prevalence_total = self.I_ac/self.population
        self.prevalence_susc = [self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]
        self.prevalence_det = [self.pI_det*self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]                         
        return
       
    def df_build(self):
        """
        Builds a dataframe with the simulation results
        """
        self.results = pd.DataFrame({'t':self.t,'dates':self.dates})
        names = ['S','E','E_d','I','I_d','R','R_d','T','T_d','Q','N']        
        aux = pd.DataFrame(np.transpose(self.sol.y),columns=names).astype(int)
        names2 = ['E_ac','I_ac','R_ac','T_ac','I_det','I_d_det','I_ac_det']
        vars2 = [self.E_ac,self.I_ac,self.R_ac,self.T_ac,self.I_det,self.I_d_det,self.I_ac_det]
        aux2 = pd.DataFrame(np.transpose(vars2),columns=names2).astype(int)

        # Prevalence
        self.prevalence = pd.DataFrame(np.transpose([self.prevalence_total,self.prevalence_susc,self.prevalence_det]),columns = ['prevalence_total','prevalence_susc','prevalence_det'])
        
        # Parameters
        nameparams = ['beta','alpha','tE_I','tI_R','rR_S','tR_S','k_Ex','k_Eacc','k_Q','k_Tr','tQ_R','tT_Q']
        beta_val = [self.beta(t) for t in self.t] 
        alpha_val = [self.alpha(t) for t in self.t] 
        tE_I_val = [self.tE_I(t) for t in self.t]
        tI_R_val = [self.tI_R(t) for t in self.t]
        rR_S_val = [self.rR_S(t) for t in self.t]
        tR_S_val = [self.tR_S(t) for t in self.t]
        k_Ex_val = [self.k_Ex(t) for t in self.t]
        k_Eacc_val = [self.k_Eacc(t) for t in self.t]
        k_Q_val = [self.k_Q(t) for t in self.t]
        k_Tr_val = [self.k_Tr(t) for t in self.t]
        tQ_R_val = [self.tQ_R(t) for t in self.t]
        tT_Q_val = [self.tT_Q(t) for t in self.t]
        self.params = pd.DataFrame(np.transpose([beta_val,alpha_val,tE_I_val,tI_R_val,rR_S_val,tR_S_val,k_Ex_val,k_Eacc_val,k_Q_val,k_Tr_val,tQ_R_val,tT_Q_val]),columns = nameparams)
        
        self.compartments = pd.concat([self.results,aux,aux2],axis=1)
        self.results = pd.concat([self.results,aux,aux2,self.params,self.prevalence],axis=1)
        self.results = self.results.astype({'S': int,'E': int,'E_d': int,'I': int,'I_d': int,'R': int,'R_d': int,'T': int,'T_d': int,'Q': int,'E_ac': int,'I_ac': int,'R_ac': int,'T_ac': int,'I_det': int,'I_d_det': int,'I_ac_det': int})

        self.resume = pd.DataFrame({'peak':int(self.peak),'peak_t':self.peak_t,'peak_date':self.peak_date},index=[0])
        return
    
    
