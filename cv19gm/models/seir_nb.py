#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIR Model
"""

import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
from datetime import timedelta

from numbalsoda import lsoda_sig, lsoda
import numba as nb


# cv19gm libraries 
import cv19gm.utils.cv19files as cv19files


class SEIR:
    """
        SEIR model object:
        Construction:
            SEIR(self, config = None, inputdata=None)

    """
    def __init__(self, config = None, inputdata=None,verbose = False, nbsolver = True, **kwargs):
        self.compartmentalmodel = "SEIR"
        
        if not config:
            #print('Missing configuration file ')
            raise('Missing configuration file')
            #return None
        
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
        
        self.nbsolver = nbsolver
        
        #if self.nbsolver:
        #    self.nbfunctions()
        
        if verbose:
            print('SEIR object created')

    # ------------------- #
    #  Valores Iniciales  #
    # ------------------- #
   
    def set_relational_values(self):
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
        
        self.N = self.popfraction*self.population
        self.S = self.N-self.E-self.I-self.R

        # External flux functions: 
        if not hasattr(self,'S_f'):
            self.S_f = lambda t:0
        
        if not hasattr(self,'E_f'):
            self.E_f = lambda t:0

        if not hasattr(self,'I_f'):
            self.I_f = lambda t:0
        
        if not hasattr(self,'R_f'):
            self.R_f = lambda t:0

        # Saturation kinetics
        if not hasattr(self,'k_I'):
            self.k_I=0

        if not hasattr(self,'k_R'):
            self.k_R=0            
    

    def set_equations(self):
        """        
        Sets Diferential Equations
        """
        # Dynamic Parameters
        alpha = nb.njit(self.alpha)
        beta = nb.njit(self.beta)
        rR_S = nb.njit(self.rR_S)
        tE_I = nb.njit(self.tE_I)
        tI_R = nb.njit(self.tI_R)
        # Flux
        S_f = nb.njit(self.S_f)
        E_f = nb.njit(self.E_f)
        I_f = nb.njit(self.I_f)
        R_f = nb.njit(self.R_f)
        
        # Static Parameters
        N = self.N
        k_I = self.k_I
        k_R = self.k_R
        
    
        # --------------------------- #
        #        Susceptibles         #
        # --------------------------- #
       
        # 0) dS/dt:
        self.dS=lambda t,S,E,I,R: S_f(t) - alpha(t)*beta(t)*S*I/(N+k_I*I + k_R*R) + rR_S(t)*R
        
        # --------------------------- #
        #           Exposed           #
        # --------------------------- #     
        
        # 1) dE/dt
        self.dE = lambda t,S,E,I,R: E_f(t) + alpha(t)*beta(t)*S*I/(N+k_I*I + k_R*R) - E/tE_I(t)
 
        # 2) Daily dE/dt
        self.dE_d = lambda t,S,E,E_d,I,R: E_f(t) + alpha(t)*beta(t)*S*I/(N+k_I*I + k_R*R) - E_d

        # --------------------------- #
        #           Infected          #
        # --------------------------- #                
        
        # 3) Active
        self.dI=lambda t,E,I: I_f(t) + E/tE_I(t) - I/tI_R(t)
        
        # 4) New Daily
        self.dI_d = lambda t,E,I_d: I_f(t) + E/tE_I(t) - I_d 

        # --------------------------- #
        #         Recovered           #
        # --------------------------- #  
        
        # 5) Total recovered
        self.dR=lambda t,I,R: R_f(t) + I/tI_R(t) - rR_S(t)*R

        # 6) Recovered per day
        self.dR_d=lambda t,I,R_d: R_f(t) + I/tI_R(t) - R_d

        # 7) External Flux:
        self.dFlux = lambda t: S_f(t) + E_f(t) + I_f(t) + R_f(t) 

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
        
        self.t = np.linspace(0, T, T+1)
        initcond = np.array([self.S,self.E,self.E_d,self.I,self.I_d,self.R,0,0]) # [S0,E0,E_d0,I0,I_d0,R0,R_d0,Flux0]
        
        if self.nbsolver:
            dS = nb.njit(self.dS)
            dE = nb.njit(self.dE)
            dE_d = nb.njit(self.dE_d)
            dI = nb.njit(self.dI)
            dI_d = nb.njit(self.dI_d)
            dR = nb.njit(self.dR)
            dR_d = nb.njit(self.dR_d)
            dFlux = nb.njit(self.dFlux)
            
            @nb.cfunc(lsoda_sig)
            def model_SEIR_graph_nb(t,y,ydot,p):
                ydot[0] = dS(t,y[0],y[1],y[3],y[5])
                ydot[1] = dE(t,y[0],y[1],y[3],y[5])
                ydot[2] = dE_d(t,y[0],y[1],y[2],y[3],y[5])
                ydot[3] = dI(t,y[1],y[3])
                ydot[4] = dI_d(t,y[1],y[4])
                ydot[5] = dR(t,y[3],y[5])
                ydot[6] = dR_d(t,y[3],y[6])
                ydot[7] = dFlux(t)
            
            funcptr = model_SEIR_graph_nb.address
            sol, success = lsoda(funcptr, initcond, self.t)
            
            self.S = sol[:,0]
            self.E = sol[:,1]
            self.E_d = sol[:,2]
            self.I = sol[:,3]
            self.I_d = sol[:,4]
            self.R = sol[:,5]
            self.R_d = sol[:,6]
            self.Flux = sol[:,7]   
            self.sol = sol         
            
        else:
            sol = solve_ivp(self.model_SEIR_graph,(t0,T), initcond,method='LSODA',t_eval=list(range(t0,T)))
            self.t=sol.t 
            self.S=sol.y[0,:]
            self.E=sol.y[1,:]
            self.E_d=sol.y[2,:]
            self.I=sol.y[3,:]
            self.I_d=sol.y[4,:]
            self.R=sol.y[5,:]
            self.R_d=sol.y[6,:]
            self.Flux=sol.y[7,:]
            
            self.sol = sol.y

        

        self.E_ac = np.cumsum(np.concatenate([np.array([self.E_ac]),self.E_d[:-1]]))
        self.I_ac = np.cumsum(np.concatenate([np.array([self.I_ac]),self.I_d[:-1]]))  
        self.R_ac = np.cumsum(np.concatenate([np.array([self.R[0]]),self.R_d[:-1]]))

        self.I_det = self.I*self.pI_det
        self.I_d_det = self.I_d*self.pI_det
        self.I_ac_det = self.I_ac*self.pI_det

        self.analytics()
        #self.df_build()
        self.solved = True
        return 
    

    def model_SEIR_graph(self,t,y):
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
       
    def df_build(self):
        """
        Builds a dataframe with the simulation results
        """
        self.results = pd.DataFrame({'t':self.t,'dates':self.dates})
        names = ['S','E','E_d','I','I_d','R','R_d','Flux']
        
        
        aux = pd.DataFrame(np.transpose(self.sol),columns=names)       

        names2 = ['E_ac','I_ac','R_ac','I_det','I_d_det','I_ac_det','prevalence_total','prevalence_susc','prevalence_det']
        vars2 = [self.E_ac,self.I_ac,self.R_ac,self.I_det,self.I_d_det,self.I_ac_det,self.prevalence_total,self.prevalence_susc,self.prevalence_det]
        aux2 = pd.DataFrame(np.transpose(vars2),columns=names2)

        self.results = pd.concat([self.results,aux,aux2],axis=1)
        self.results = self.results.astype({'S': int,'E': int,'E_d': int,'I': int,'I_d': int,'R': int,'R_d': int,'E_ac': int,'I_ac': int,'R_ac': int,'I_det': int,'I_d_det': int,'I_ac_det': int})

        self.resume = pd.DataFrame({'peak':int(self.peak),'peak_t':self.peak_t,'peak_date':self.peak_date},index=[0])
        return

    
    def nbfunctions(self):
        print('Using numbaLSODA')
        # Equations
        #self.dS = nb.njit(self.dS)
        #self.dE =  nb.njit(self.dE)
        #self.dE_d = nb.njit(self.dE_d)
        #self.dI =  nb.njit(self.dI)
        #self.dI_d =  nb.njit(self.dI_d)
        #self.dR =  nb.njit(self.dR)
        #self.dR_d =  nb.njit(self.dR_d)
        #self.dFlux =  nb.njit(self.dFlux)
        # Parameters
        self.alpha = nb.njit(self.alpha)
        self.beta = nb.njit(self.beta)
        self.rR_S = nb.njit(self.rR_S)
        self.tE_I = nb.njit(self.tE_I)
        self.tI_R = nb.njit(self.tI_R)
        # Flux
        self.S_f = nb.njit(self.S_f)
        self.E_f = nb.njit(self.E_f)
        self.I_f = nb.njit(self.I_f)
        self.R_f = nb.njit(self.R_f)
        
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

        
