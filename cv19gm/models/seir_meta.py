#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIR Meta-populations Model
"""

import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd

import cv19gm.utils.cv19files as cv19files
import cv19gm.utils.cv19mobility as cv19mobility

""" 
ToDo
    * Add accumulated variables
    * Simplify results_build using params_build 
    * Optimize code (remove unnecessary variables)

"""

class SEIRMETA:  
    """|
        SEIRMETA model object:
        Construction:
            SEIRMETA(self, config = None)

    """
    def __init__(self, config = None, verbose = False, Phi = None, Phi_T = None, seed=None, method = 0, **kwargs):    
        if verbose and not config:
            print('Warning: Using default configuration file')
         
        self.compartmentalmodel = "SEIR_Metapopulation"
        self.kwargs = kwargs        
        self.method = method
        # ------------------------------- #
        #         Parameters Load         #
        # ------------------------------- #
        self.config = config
        if verbose:
            print('Loading configuration file')
                      
        cv19files.loadconfig(self,config,**kwargs)
        
        # Mobility matrix
        if Phi:
            self.Phi = Phi
            if Phi_T:
                self.Phi_T = Phi_T
            else:
                self.Phi_T = lambda t: Phi(t).transpose()
                
        else:
            if verbose:
                print('Warning: Missing human mobility matrix, using a random matrix instead')
            self.Phi, self.Phi_T = cv19mobility.create_dynamic_mobility(mobility_model='random', dynamic_pattern='symmetric',populations=self.population, seed=seed, transposed = True)
            
        if verbose:
            print('Initializing parameters and variables')
        self.set_initial_values()
        if verbose:
            print('Building equations')          
        self.set_equations()

        self.solved = False
        if verbose:
            print('SEIR object created')

    def set_initial_values(self):
        # Exposed
        #The np.array cast is transitory hopefuly:
        if not hasattr(self,'E') or not self.E:
            self.E = np.array(self.mu)*np.array(self.I)
            self.E_d = np.array(self.mu)*np.array(self.I_d)

        if not hasattr(self,'E_ac') or not self.E_ac:
            self.E_ac= 0 
       
        # Valores globales
        if not hasattr(self,'popfraction') or not self.popfraction:
            self.popfraction = 1
        
        self.nregions = len(self.population)
        
        self.N = self.popfraction*np.array(self.population)
        self.S = self.N-self.E-self.I-self.R
        self.nodes = len(self.population) # Amount of nodes /meta-populations
    
    def set_equations(self):
        """        
        Sets Diferential Equations
        """
        # --------------------------- #
        #        Susceptibles         #
        # --------------------------- #
       
        # 0) dS/dt:
        self.dS=lambda t,S,I,R,N: - np.dot(np.diag(S*I/N),(self.alpha(t)*self.beta(t))) + self.rR_S(t)*R  + self.phi_S(t,S,N)
        
        # --------------------------- #
        #           Exposed           #
        # --------------------------- #     
        
        # 1) dE/dt
        self.dE = lambda t,S,E,I,N: np.dot(np.diag(S*I/N),(self.alpha(t)*self.beta(t))) - E/self.tE_I(t) + self.phi_E(t,E,N)
 
        # 2) Daily dE/dt
        self.dE_d = lambda t,S,E,E_d,I,N: np.dot(np.diag(S*I/N),(self.alpha(t)*self.beta(t))) + self.phi_E(t,E,N) - E_d

        # --------------------------- #
        #           Infected          #
        # --------------------------- #                
        
        # 3) Active
        self.dI=lambda t,E,I,N: E/self.tE_I(t) - I/self.tI_R(t) + self.phi_I(t,I,N)
        
        # 4) New Daily
        self.dI_d = lambda t,E,I,I_d,N: E/self.tE_I(t) + self.phi_I(t,I,N) - I_d 

        # --------------------------- #
        #         Recovered           #
        # --------------------------- #  
        
        # 5) Total recovered
        self.dR=lambda t,I,R,N: I/self.tI_R(t) - self.rR_S(t)*R + self.phi_R(t,R,N)

        # 6) Recovered per day
        self.dR_d=lambda t,I,R,R_d,N: I/self.tI_R(t) + self.phi_R(t,R,N) - R_d

        # 7) Population Flux:
        self.dN = lambda t,S,E,I,R,N: self.phi_S(t,S,N) + self.phi_E(t,E,N) + self.phi_I(t,I,N) + self.phi_R(t,R,N)
        
        # --------------------------- #
        #         People Flux         #
        # --------------------------- # 
        # Method 0
        # Original system of equations
        if self.method == 0:
            np_ones = np.ones(self.nregions)
            self.phi_S = lambda t,S,N: np.dot(self.Phi_T(t),(S/N)) - np.dot(np.dot(np.diag(S/N),self.Phi(t)),np_ones)
            self.phi_E = lambda t,E,N: np.dot(self.Phi_T(t),(E/N)) - np.dot(np.dot(np.diag(E/N),self.Phi(t)),np_ones)
            self.phi_I = lambda t,I,N: np.dot(self.Phi_T(t),(I/N)) - np.dot(np.dot(np.diag(I/N),self.Phi(t)),np_ones)
            self.phi_R = lambda t,R,N: np.dot(self.Phi_T(t),(R/N)) - np.dot(np.dot(np.diag(R/N),self.Phi(t)),np_ones)       
       
                
	    # Method 1 
        # This method replaces Phi_t = Phi(t), which is not correct, but it's much faster and I don't understand why
        elif self.method == 1:
            np_ones = np.ones(self.nregions)
            self.phi_S = lambda t,S,N: np.dot(self.Phi(t),(S/N)) - np.dot(np.dot(np.diag(S/N),self.Phi(t)),np_ones)
            self.phi_E = lambda t,E,N: np.dot(self.Phi(t),(E/N)) - np.dot(np.dot(np.diag(E/N),self.Phi(t)),np_ones)
            self.phi_I = lambda t,I,N: np.dot(self.Phi(t),(I/N)) - np.dot(np.dot(np.diag(I/N),self.Phi(t)),np_ones)
            self.phi_R = lambda t,R,N: np.dot(self.Phi(t),(R/N)) - np.dot(np.dot(np.diag(R/N),self.Phi(t)),np_ones)
    
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
        self.R_d = np.zeros(self.nregions)
        initcond = np.array([self.S,self.E,self.E_d,self.I,self.I_d,self.R,self.R_d,self.N]).flatten()
        
        sol = solve_ivp(self.solver_equations,(t0,T), initcond,method=method,t_eval=list(range(t0,T)))
        
        self.sol = sol
        self.t=sol.t         
        sol.y = sol.y.reshape([8,self.nregions,len(self.t)])

        self.S=sol.y[0]
        self.E=sol.y[1]
        self.E_d=sol.y[2]
        self.I=sol.y[3]
        self.I_d=sol.y[4]
        self.R=sol.y[5]
        self.R_d=sol.y[6]
        self.N=sol.y[7]

        #self.E_ac = np.cumsum(self.E_d)
        #self.I_ac = np.cumsum(self.I_d) + self.I_ac # second term is the initial condition
        #self.R_ac = np.cumsum(self.R_d)

        self.results_build()
        self.global_results_build()
        self.solved = True

        return 

    def solver_equations(self,t,y):
        y = y.reshape(8, self.nregions) #8 is the number of equations
        ydot=np.zeros(np.shape(y))
        ydot[0]=self.dS(t,y[0],y[3],y[5],y[7])
        ydot[1]=self.dE(t,y[0],y[1],y[3],y[7])
        ydot[2]=self.dE_d(t,y[0],y[1],y[2],y[3],y[7])        
        ydot[3]=self.dI(t,y[1],y[3],y[7])
        ydot[4]=self.dI_d(t,y[1],y[3],y[4],y[7])        
        ydot[5]=self.dR(t,y[3],y[5],y[7])
        ydot[6]=self.dR_d(t,y[3],y[5],y[6],y[7])        
        ydot[7]=self.dN(t,y[0],y[1],y[3],y[5],y[7])                                        
        return(ydot.flatten())

    def results_build(self):
        """
        Params shouldn't be int! 
        Builds a dataframe with the simulation results and parameters 
        Output structure: 
        't','S','E','E_d','I','I_d','R','R_d','beta','tE_I','tI_R','rR_S','node'
         0, ...
         1, ...
          
        """
        names = ['t','S','E','E_d','I','I_d','R','R_d','alpha','beta','tE_I','tI_R','rR_S','node']
        
        # Parameters
        alpha_val = [[self.alpha(t)[j] for t in self.t] for j in range(self.nodes)]
        beta_val = [[self.beta(t)[j] for t in self.t] for j in range(self.nodes)]        
        tE_I_val = [self.tE_I(t) for t in self.t]
        tI_R_val = [self.tI_R(t) for t in self.t]
        rR_S_val = [self.rR_S(t) for t in self.t]
        
        self.results = []
        for i in range(self.nodes):
            node = [i]*len(self.t)
            self.results.append(pd.DataFrame(dict(zip(names,[self.t,self.S[i],self.E[i],self.E_d[i],self.I[i],self.I_d[i],self.R[i],self.R_d[i],alpha_val[i],beta_val[i],tE_I_val,tI_R_val,rR_S_val,node]))))        
        self.results = pd.concat(self.results,ignore_index=True).astype(int)
        return
    
    def global_results_build(self):
        """Agregated results data frame
        """
        self.S_tot = self.S.sum(axis=0)
        self.E_tot = self.E.sum(axis=0)
        self.E_d_tot = self.E_d.sum(axis=0)
        self.I_tot = self.I.sum(axis=0)
        self.I_d_tot = self.I_d.sum(axis=0)
        self.R_tot = self.R.sum(axis=0)
        self.R_d_tot = self.R_d.sum(axis=0)
        
        names = ['t','S','E','E_d','I','I_d','R','R_d']        
        self.global_results = pd.DataFrame(np.array([self.t,self.S_tot,self.E_tot,self.E_d_tot,self.I_tot,self.I_d_tot,self.R_tot,self.R_d_tot]).transpose(),columns=names).astype(int)        
        return        
        
    def params_df_build(self):
        """
        Builds a dataframe with the simulation parameters over time
        """
        
        names = ['t','alpha','beta','tE_I','tI_R','rR_S','node']  
        alpha_val = [[self.alpha(t)[j] for t in self.t] for j in range(self.nodes)]               
        beta_val = [[self.beta(t)[j] for t in self.t] for j in range(self.nodes)]        
        tE_I_val = [self.tE_I(t) for t in self.t]
        tI_R_val = [self.tI_R(t) for t in self.t]
        rR_S_val = [self.rR_S(t) for t in self.t]

        self.params = []
        for i in range(self.nodes):
            node = [i]*len(self.t)
            self.params.append(pd.DataFrame(dict(zip(names,[self.t,alpha_val[i],beta_val[i],tE_I_val,tI_R_val,rR_S_val,node]))))
        
        self.params = pd.concat(self.params,ignore_index=True)
        return