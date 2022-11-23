#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIR Meta-population Model
"""

from unittest import result
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
from datetime import timedelta

# cv19gm libraries 
#import os
#import sys
#path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
#sys.path.insert(1, path)

#import data.cv19data as cv19data
#import utils.cv19timeutils as cv19timeutils
#import utils.cv19functions as cv19functions
import cv19gm.utils.cv19files as cv19files
import cv19gm.utils.cv19mobility as cv19mobility

""" To Do
* Y tener precargaada la matriz de movilidad como transpuesta
* Matriz de movilidad de función a tensor

"""

class SEIRMETA:  
    """|
        SEIR model object:
        Construction:
            SEIR(self, config = None, inputdata=None)

    """
    def __init__(self, config = None, inputdata=None,verbose = False, Phi = None, Phi_T = None, seed=None, method = 0, **kwargs):    
        if not config:
            #print('Missing configuration file ')
            raise('Missing configuration file')
            #return None
        self.kwargs = kwargs
        
        self.method = method
        # ------------------------------- #
        #         Parameters Load         #
        # ------------------------------- #
        self.config = config
        if verbose:
            print('Loading configuration file')          
        cv19files.loadconfig(self,config,inputdata,**kwargs)
        
        # Definition of mobility matrix (Work in progress, it will be done inside the cb19mobility lib)
        if Phi:
            self.Phi = Phi
            if Phi_T:
                self.Phi_T = Phi_T
            else:
                self.Phi_T = lambda t: Phi(t).transpose()
                
        else:
            print('Missing flux dynamics, using a random matrix instead')
            self.Phi, self.Phi_T = cv19mobility.rnd_flux_symmetric(self.population,seed=seed, transposed = True)
            
        #if not hasattr(self,'Phi') or not self.Phi:
        #    print('Missing flux dynamics, using a random matrix instead')
        #    self.Phi = cv19mobility.rnd_flux_symmetric(self.population)
            
        if verbose:
            print('Initializing parameters and variables')
        self.set_initial_values()
        if verbose:
            print('Building equations')          
        self.set_equations()

        self.solved = False
        if verbose:
            print('SEIR object created')

    # ------------------- #
    #  Valores Iniciales  #
    # ------------------- #
   
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
        self.dS=lambda t,S,I,R,N: - np.diag(S*I/N)@self.beta(t) + self.rR_S(t)*R  + self.phi_S(t,S,N)
        
        # --------------------------- #
        #           Exposed           #
        # --------------------------- #     
        
        # 1) dE/dt
        self.dE = lambda t,S,E,I,N: np.diag(S*I/N)@self.beta(t) - E/self.tE_I(t) + self.phi_E(t,E,N)
 
        # 2) Daily dE/dt
        self.dE_d = lambda t,S,E,E_d,I,N: np.diag(S*I/N)@self.beta(t) + self.phi_E(t,E,N) - E_d

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
        
        if self.method == 0:        
        # Version original
            self.phi_S = lambda t,S,N: self.Phi(t).transpose()@(S/N) - np.diag(S/N)@self.Phi(t)@np.ones(self.nregions)
            self.phi_E = lambda t,E,N: self.Phi(t).transpose()@(E/N) - np.diag(E/N)@self.Phi(t)@np.ones(self.nregions)
            self.phi_I = lambda t,I,N: self.Phi(t).transpose()@(I/N) - np.diag(I/N)@self.Phi(t)@np.ones(self.nregions)
            self.phi_R = lambda t,R,N: self.Phi(t).transpose()@(R/N) - np.diag(R/N)@self.Phi(t)@np.ones(self.nregions)

        elif self.method == 1:        
        # Primera propuesta
            np_ones = np.ones(self.nregions)
            self.phi_S = lambda t,S,N: self.Phi(t).transpose()@(S/N) - np.diag(S/N)@self.Phi(t)@np_ones
            self.phi_E = lambda t,E,N: self.Phi(t).transpose()@(E/N) - np.diag(E/N)@self.Phi(t)@np_ones
            self.phi_I = lambda t,I,N: self.Phi(t).transpose()@(I/N) - np.diag(I/N)@self.Phi(t)@np_ones
            self.phi_R = lambda t,R,N: self.Phi(t).transpose()@(R/N) - np.diag(R/N)@self.Phi(t)@np_ones
        
        # segunda propuesta 
        elif self.method == 2:
            self.Phi_matrix = cv19mobility.mobility_to_tensor(self.Phi,self.tsim)
            self.Phi_matrix_T = cv19mobility.mobility_transposed(self.Phi_matrix)
            np_ones = np.ones(self.nregions)
        
            self.phi_S = lambda t,S,N: self.Phi_matrix_T[int(2*t)]@(S/N) - np.diag(S/N)@self.Phi_matrix[int(2*t)]@np_ones
            self.phi_E = lambda t,E,N: self.Phi_matrix_T[int(2*t)]@(E/N) - np.diag(E/N)@self.Phi_matrix[int(2*t)]@np_ones
            self.phi_I = lambda t,I,N: self.Phi_matrix_T[int(2*t)]@(I/N) - np.diag(I/N)@self.Phi_matrix[int(2*t)]@np_ones
            self.phi_R = lambda t,R,N: self.Phi_matrix_T[int(2*t)]@(R/N) - np.diag(R/N)@self.Phi_matrix[int(2*t)]@np_ones        
        
        # Tercera propuesta 
        elif self.method == 3:
            self.Phi_matrix = cv19mobility.mobility_to_tensor(self.Phi,self.tsim)
            self.Phi_matrix_T = cv19mobility.mobility_transposed(self.Phi_matrix)
            np_ones = np.ones(self.nregions)

            self.phi_S = lambda t,S,N: np.dot(self.Phi_matrix_T[int(2*t)],(S/N)) - np.dot(np.dot(np.diag(S/N),self.Phi_matrix[int(2*t)]),np_ones)
            self.phi_E = lambda t,E,N: np.dot(self.Phi_matrix_T[int(2*t)],(E/N)) - np.dot(np.dot(np.diag(E/N),self.Phi_matrix[int(2*t)]),np_ones)
            self.phi_I = lambda t,I,N: np.dot(self.Phi_matrix_T[int(2*t)],(I/N)) - np.dot(np.dot(np.diag(I/N),self.Phi_matrix[int(2*t)]),np_ones)
            self.phi_R = lambda t,R,N: np.dot(self.Phi_matrix_T[int(2*t)],(R/N)) - np.dot(np.dot(np.diag(R/N),self.Phi_matrix[int(2*t)]),np_ones)
                
        # Tercera propuesta 
        elif self.method == 4:
            np_ones = np.ones(self.nregions)
            self.phi_S = lambda t,S,N: np.dot(self.Phi_T(t),(S/N)) - np.dot(np.dot(np.diag(S/N),self.Phi(t)),np_ones)
            self.phi_E = lambda t,E,N: np.dot(self.Phi_T(t),(E/N)) - np.dot(np.dot(np.diag(E/N),self.Phi(t)),np_ones)
            self.phi_I = lambda t,I,N: np.dot(self.Phi_T(t),(I/N)) - np.dot(np.dot(np.diag(I/N),self.Phi(t)),np_ones)
            self.phi_R = lambda t,R,N: np.dot(self.Phi_T(t),(R/N)) - np.dot(np.dot(np.diag(R/N),self.Phi(t)),np_ones)

                    
    def integrate(self,t0=0,T=None,h=0.01):
        print('The use of integrate() is now deprecated. Use solve() instead.')
        self.solve(t0=t0,T=T,h=h)

    def run(self,t0=0,T=None,h=0.01):
        #print('The use of integrate() is now deprecated. Use solve() instead.')
        self.solve(t0=t0,T=T,h=h)

    # Scipy
    def solve(self,t0=0,T=None,h=0.01):
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
        initcond = np.array([self.S,self.E,self.E_d,self.I,self.I_d,self.R,self.R_d,self.N]).flatten() # [S0,E0,E_d0,I0,I_d0,R0,R_d0,Flux0]
        
        sol = solve_ivp(self.model_SEIR_graph,(t0,T), initcond,method='LSODA',t_eval=list(range(t0,T)))
        
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

        #self.analytics()
        self.results_build()
        self.global_results_build()
        #self.underreport()
        self.solved = True

        return 

    def model_SEIR_graph(self,t,y):
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
        Builds a dataframe with the simulation results and parameters 
        Output structure: 
        't','S','E','E_d','I','I_d','R','R_d','beta','tE_I','tI_R','rR_S','node'
         0, ...
         1, ...
          
        """
        names = ['t','S','E','E_d','I','I_d','R','R_d','beta','tE_I','tI_R','rR_S','node']
        # Parameters
        
        beta_val = [[self.beta(t)[j] for t in self.t] for j in range(self.nodes)]        
        tE_I_val = [self.tE_I(t) for t in self.t]
        tI_R_val = [self.tI_R(t) for t in self.t]
        rR_S_val = [self.rR_S(t) for t in self.t]
        
        self.results = []
        for i in range(self.nodes):
            node = [i]*len(self.t)
            self.results.append(pd.DataFrame(dict(zip(names,[self.t,self.S[i],self.E[i],self.E_d[i],self.I[i],self.I_d[i],self.R[i],self.R_d[i],beta_val[i],tE_I_val,tI_R_val,rR_S_val,node]))))        
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
        
        names = ['t','beta','tE_I','tI_R','rR_S','node']                 
        beta_val = [[self.beta(t)[j] for t in self.t] for j in range(self.nodes)]
        
        tE_I_val = [self.tE_I(t) for t in self.t]
        tI_R_val = [self.tI_R(t) for t in self.t]
        rR_S_val = [self.rR_S(t) for t in self.t]

        self.params = []
        for i in range(self.nodes):
            node = [i]*len(self.t)
            self.params.append(pd.DataFrame(dict(zip(names,[self.t,beta_val[i],tE_I_val,tI_R_val,rR_S_val,node]))))
        
        self.params = pd.concat(self.params,ignore_index=True)
        return
    
    

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
            self.peak_date = self.initdate+timedelta(days=int(round(self.peak_t)))
        else:
            self.dates = [None for i in range(len(self.t))]
            self.peak_date = None            

        # Prevalence: 
        self.prevalence_total = self.I_ac/self.population
        self.prevalence_susc = [self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]
        self.prevalence_det = [self.pI_det*self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))]                         
        return
           

    """ 
    def calculateindicators(self):
        self.R_ef
        self.SHFR
        # SeroPrevalence Calculation
        # Errors (if real data)
        # Active infected
        print('wip')

    def resume(self):        
        print("Resumen de resultados:")
        qtype = ""
        for i in range(self.numescenarios):
            if self.inputarray[i][-1]==0:
                qtype = "Cuarentena total"
            if self.inputarray[i][-1]>0:
                qtype ="Cuarentena Dinámica"            

            print("Escenario "+str(i))
            print("Tipo de Cuarentena: "+qtype+'\nmov_rem: '+str(self.inputarray[i][2])+'\nmov_max: '+str(self.inputarray[i][2])+
            "\nInicio cuarentena: "+(self.initdate+timedelta(days=self.inputarray[i][4])).strftime('%Y/%m/%d')+"\nFin cuarentena: "+(self.initdate+timedelta(days=self.inputarray[i][5])).strftime('%Y/%m/%d'))
            print("Peak infetados \n"+"Peak value: "+str(self.peak[i])+"\nPeak date: "+str(self.peak_date[i]))
            print("Fallecidos totales:"+str(max(self.B[i])))
            print("Fecha de colapso hospitalario \n"+"Camas: "+self.H_colapsedate[i]+"\nVentiladores: "+self.V_colapsedate[i])
            print("\n")
    """

        
