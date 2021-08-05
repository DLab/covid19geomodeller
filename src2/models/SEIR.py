#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIR Model
"""

import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
import toml
from datetime import datetime
from datetime import timedelta


# Deprecated
#from scipy.special import expit
#from joblib import Parallel, delayed
#import multiprocessing  
#from scipy import signal
#from numpy import linalg as LA 

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
    def __init__(self, config = None, inputdata=None,aux=None,**kwargs):
        
        if aux:
            kwargs.update(aux)
        # Parameter that defines sigmoid's raising speed
        self.gw=20
        self.config = config
        
        if not config:
            print('Missing configuration file ')
            return None

        # ------------------------------- #
        #         Parameters Load         #
        # ------------------------------- #
        cv19files.loadconfig(self,config,inputdata,**kwargs)
        self.setrelationalvalues()
        self.setequations()

        self.solved = False
        #print('SEIR object created')

    # ------------------- #
    #  Valores Iniciales  #
    # ------------------- #
   
    def setrelationalvalues(self):
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
                    

    def setequations(self):
        """
        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #
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



    def integr_sci(self,t0=0,T=None,h=0.01):
        self.integrate(t0=0,T=None,h=0.01)
        return

    # Scipy
    def integrate(self,t0=0,T=None,h=0.01):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,mu*(self.I_ac)T] function doesn't
        #start. Or it's start if class object is initialze.

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
            
        elif False:#((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            E0 = self.E[idx]
            E_d0 = self.E_d[idx]
            E_ac0 = self.E_ac[idx]

            S0=self.S[idx]
            Ias0=self.Ias[idx]
            Imi0=self.Imi[idx]
            Ise0=self.Ise[idx]
            Icr0=self.Icr[idx]

            Ias_d0=self.Ias_d[idx]
            Imi_d0=self.Imi_d[idx]
            Ise_d0=self.Ise_d[idx]
            Icr_d0=self.Icr_d[idx]

            Ias_ac0=self.Ias_ac[idx]
            Imi_ac0=self.Imi_ac[idx]
            Ise_ac0=self.Ise_ac[idx]
            Icr_ac0=self.Icr_ac[idx]
        
            Hse0=self.Hse[idx]
            Hout0=self.Hout[idx]
            V0=self.V[idx]

            Hse_d0=self.Hse_d[idx]
            Hout_d0=self.Hout_d[idx]
            V_d0=self.V_d[idx]

            Hse_ac0=self.Hse_ac[idx]
            Hout_ac0=self.Hout_ac[idx]
            V_ac0=self.V_ac[idx]

            R0=self.R[idx]
            R_d0=self.R_d[idx]

            D0=self.D[idx]
            B0=self.B[idx]

            Ise_D_d0 = self.Ise_D_d[idx]
            Icr_D_d0 = self.Icr_D_d[idx]
            Hse_D_d0 = self.Hse_D_d[idx]
            V_D_d0 = self.V_D_d[idx]

            Ise_D_ac0 = self.Ise_D_ac[idx]
            Icr_D_ac0 = self.Icr_D_ac[idx]
            Hse_D_ac0 = self.Hse_D_ac[idx]
            V_D_ac0 = self.V_D_ac[idx] 

            V_need0=self.V[idx]
            Hse_need0=self.Hse[idx]
            Hout_need0= self.Hout[idx]

            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)

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
        self.dfbuild(sol)
        self.solved = True

        return(sol)





    # sckits: slower but better
    def integr(self,t0=0,T=None,h=0.01,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.
        print('Import scikits-odes')
        from scikits.odes.odeint import odeint

        if T is None:
            T = self.tsim

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
            
        elif False:#((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            S0=self.S[idx]
            E0=self.E[idx]
            
            I0=self.I[idx]
            R0=self.R[idx]                        
            I_ac0=self.I_ac[idx]
            I_d0=self.I_d[idx]

            e0 = self.e[idx]
            e_I0 = self.e_I[idx]
            
            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)

        else:
            print('Already solved')
            return()

        
        def model_SEIR_graph(t,y,ydot):
            
            ydot[0]=self.dS(t,y[0],y[1],y[3],y[5])

            ydot[1]=self.dE(t,y[0],y[1],y[3],y[5])
            ydot[2]=self.dE_d(t,y[0],y[1],y[2],y[3],y[5])

            ydot[3]=self.dI(t,y[1],y[3])
            ydot[4]=self.dI_d(t,y[1],y[4])

            ydot[5]=self.dR(t,y[3],y[5])
            ydot[6]=self.dR_d(t,y[3],y[6])

            ydot[7]=self.dFlux(t)

            
        initcond = np.array([S0,E0,E_d0,I0,I_d0,R0,R_d0,Flux0])                                


        sol = odeint(model_SEIR_graph, self.t, initcond,method='admo')
        
        self.sol = sol
        self.t=sol.values.t 
        
        self.S=sol.values.y[:,0]
        self.E=sol.values.y[:,1]        
        self.E_d=sol.values.y[:,2]
        self.I=sol.values.y[:,3]
        self.I_d=sol.values.y[:,4]
        self.R=sol.values.y[:,5]
        self.R_d=sol.values.y[:,6]
        self.Flux=sol.values.y[:,7]

        self.E_ac = np.cumsum(self.E_d)
        self.I_ac = np.cumsum(self.I_d) + self.I_ac
        self.R_ac = np.cumsum(self.R_d)

        self.I_det = self.I*self.pI_det
        self.I_d_det = self.I_d*self.pI_det
        self.I_ac_det = self.I_ac*self.pI_det

        self.analytics()
        self.dfbuild(sol)

        self.solved = True
        
        return(sol)

    def analytics(self):
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
       
    def dfbuild(self,sol):        
        self.results = pd.DataFrame({'t':self.t,'dates':self.dates})
        names = ['S','E','E_d','I','I_d','R','R_d','Flux']
        
        self.aux = pd.DataFrame(np.transpose(sol.y),columns=names)       

        names2 = ['E_ac','I_ac','R_ac','I_det','I_d_det','I_ac_det','prevalence_total','prevalence_susc','prevalence_det']
        vars2 = [self.E_ac,self.I_ac,self.R_ac,self.I_det,self.I_d_det,self.I_ac_det,self.prevalence_total,self.prevalence_susc,self.prevalence_det]
        self.aux2 = pd.DataFrame(np.transpose(vars2),columns=names2)

        self.results = pd.concat([self.results,self.aux,self.aux2],axis=1)
        self.results = self.results.astype({'S': int,'E': int,'E_d': int,'I': int,'I_d': int,'R': int,'R_d': int,'E_ac': int,'I_ac': int,'R_ac': int,'I_det': int,'I_d_det': int,'I_ac_det': int})

        self.resume = pd.DataFrame({'peak':int(self.peak),'peak_t':self.peak_t,'peak_date':self.peak_date},index=[0])



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

        