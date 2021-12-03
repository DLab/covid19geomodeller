#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIR Model
"""

import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
from datetime import datetime
from datetime import timedelta

# cv19gm libraries 
import os
import sys
path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(1, path)

import data.cv19data as cv19data
import utils.cv19timeutils as cv19timeutils
import utils.cv19functions as cv19functions
import utils.cv19files as cv19files



class SEIR:  
    """
        SEIR model object:
        Construction:
            SEIR(self, config = None, inputdata=None)

    """
    def __init__(self, config = None, inputdata=None,verbose = False, **kwargs):
    
        if not config:
            print('Missing configuration file ')
            return None
        
        # ------------------------------- #
        #         Parameters Load         #
        # ------------------------------- #
        self.config = config
        if verbose:
            print('Loading configuration file')          
        cv19files.loadconfig(self,config,inputdata,**kwargs)
        if verbose:
            print('Initializing parameters and variables')
        self.set_relational_values()
        if verbose:
            print('Building equations')          
        self.set_equations()

        self.solved = False
        if verbose:
            print('SEIR object created')

    # ------------------- #
    #  Valores Iniciales  #
    # ------------------- #
   
    def set_relational_values(self):
        # Active infected
        if self.I_det:
            self.I = self.I_det/self.pI_det
        else:
            self.I_det = self.I*self.pI_det


        # New daily Infected
        if self.I_d_det:
            self.I_d = self.I_d_det/self.pI_det
        else:
            self.I_d_det = self.I_d*self.pI_det


        # Accumulated Infected
        if self.I_ac_det:
            self.I_ac = self.I_ac_det/self.pI_det
        else:
            self.I_ac_det = self.I_ac*self.pI_det

        
        # Exposed
        #if not self.Einit:
        self.E = self.mu*self.I
        self.E_d=self.mu*self.I_d                
        self.E_ac=self.mu*self.I_ac
       
        # Valores globales
        self.N = self.seroprevfactor*self.population
        self.S = self.N-self.E-self.I-self.R
                    

    def set_equations(self):
        """        
        Sets Diferential Equations
        """
        # --------------------------- #
        #        Susceptibles         #
        # --------------------------- # 
       
        # 0) dS/dt:
        self.dS=lambda t,S,E,I,R: self.S_f(t) - self.alpha(t)*self.beta(t)*S*(self.expinfection*E+I)/(self.N+self.k_I*I + self.k_R*R) + self.rR_S(t)*R
        
        # --------------------------- #
        #           Exposed           #
        # --------------------------- #        
        
        # 1) dE/dt
        self.dE = lambda t,S,E,I,R: self.E_f(t) + self.alpha(t)*self.beta(t)*S*(self.expinfection*E+I)/(self.N+self.k_I*I + self.k_R*R) - E/self.tE_I(t)
 
        # 2) Daily dE/dt
        self.dE_d = lambda t,S,E,E_d,I,R: self.E_f(t) + self.alpha(t)*self.beta(t)*S*(self.expinfection*E+I)/(self.N+self.k_I*I + self.k_R*R) - E_d
        
        # Accumulated dE/dt
        #self.dE_ac = lambda t,S,E,I,R: self.alpha(t)*self.beta(t)*S*(self.expinfection*E+I)/(self.N+self.k_I*I + self.k_R*R)


        # --------------------------- #
        #           Infected          #
        # --------------------------- #                
        
        # 3) Active
        self.dI=lambda t,E,I: self.I_f(t) + E/self.tE_I(t) - I/self.tI_R(t)
        
        # 4) New Daily
        self.dI_d = lambda t,E,I_d: self.I_f(t) + E/self.tE_I(t) - I_d 

        # Accummulated
        #self.dI_ac = lambda t,E: self.pE_I/self.tE_I*E

        # --------------------------- #
        #         Recovered           #
        # --------------------------- #  
        
        # 5) Total recovered
        self.dR=lambda t,I,R: self.R_f(t) + I/self.tI_R(t) - self.rR_S(t)*R

        # 6) Recovered per day
        self.dR_d=lambda t,I,R_d: self.R_f(t) + I/self.tI_R(t) - R_d

        # Recovered Accumulated
        #self.dR_ac =lambda t,I: self.pI_R/self.tI_R(t)*I 

        # 7) External Flux:
        self.dFlux = lambda t: self.S_f(t) + self.E_f(t) + self.I_f(t) + self.R_f(t) 

    def integrate(self,t0=0,T=None,h=0.01):
        print('The use of integrate() is now deprecated. Use solve() instead.')
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
        if not self.solved:

            S0=self.S
            if self.E:
                E0 = self.E
                E_d0 = self.E_d
            else:
                E0 = self.mu*(self.I)
                E_d0 = self.mu*(self.I_d)
            I0=self.I
            I_d0=self.I_d
            R0=self.R
            R_d0=0

            Flux0=0            

            self.t=np.arange(t0,T+h,h)
            

        else:
            #print('Already solved')
            return()
            

        initcond = np.array([S0,E0,E_d0,I0,I_d0,R0,R_d0,Flux0])

        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dS(t,y[0],y[1],y[3],y[5])

            ydot[1]=self.dE(t,y[0],y[1],y[3],y[5])
            ydot[2]=self.dE_d(t,y[0],y[1],y[2],y[3],y[5])

            ydot[3]=self.dI(t,y[1],y[3])
            ydot[4]=self.dI_d(t,y[1],y[4])

            ydot[5]=self.dR(t,y[3],y[5])
            ydot[6]=self.dR_d(t,y[3],y[6])

            ydot[7]=self.dFlux(t)      
                                          
            return(ydot)

         
        
        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA',t_eval=list(range(t0,T)))
        
        self.sol = sol
        self.t=sol.t 
        
        self.S=sol.y[0,:]
        self.E=sol.y[1,:]
        self.E_d=sol.y[2,:]
        self.I=sol.y[3,:]
        self.I_d=sol.y[4,:]
        self.R=sol.y[5,:]
        self.R_d=sol.y[6,:]
        self.Flux=sol.y[7,:]

        self.E_ac = np.cumsum(self.E_d)
        self.I_ac = np.cumsum(self.I_d) + self.I_ac
        self.R_ac = np.cumsum(self.R_d)

        self.I_det = self.I*self.pI_det
        self.I_d_det = self.I_d*self.pI_det
        self.I_ac_det = self.I_ac*self.pI_det

        self.analytics()
        self.df_build()
        self.solved = True

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
        names = ['S','E','E_d','I','I_d','R','R_d','Flux']
        
        aux = pd.DataFrame(np.transpose(self.sol.y),columns=names)       

        names2 = ['E_ac','I_ac','R_ac','I_det','I_d_det','I_ac_det','prevalence_total','prevalence_susc','prevalence_det']
        vars2 = [self.E_ac,self.I_ac,self.R_ac,self.I_det,self.I_d_det,self.I_ac_det,self.prevalence_total,self.prevalence_susc,self.prevalence_det]
        aux2 = pd.DataFrame(np.transpose(vars2),columns=names2)

        self.results = pd.concat([self.results,aux,aux2],axis=1)
        self.results = self.results.astype({'S': int,'E': int,'E_d': int,'I': int,'I_d': int,'R': int,'R_d': int,'E_ac': int,'I_ac': int,'R_ac': int,'I_det': int,'I_d_det': int,'I_ac_det': int})

        self.resume = pd.DataFrame({'peak':int(self.peak),'peak_t':self.peak_t,'peak_date':self.peak_date},index=[0])
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

        
