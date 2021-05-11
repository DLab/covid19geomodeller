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
from datetime import datetime
from datetime import timedelta
import toml

# cv19gm libraries 
import sys
from pathlib import Path
#sys.path.insert(1, '/utils/')
sys.path.insert(1, '../utils/')

import cv19data
import cv19timeutils
import cv19functions

"""
To do:
  - Ejecutar setparams sólo si es que no se usa archivo de configuración
  - Estandarizar nombres de variables (creo que está listo)
  - Función para construir condiciones iniciales y asi optimizar el código
  - Create reports function inside class

"""


class SEIR:  
    """
        SEIRHVD Object:
        Construction:
            SEIRHVD(tsim,beta,mu,alpha,k=0,H_cap=30,V_cap=20)
        with:
            tsim: Simulation time
            beta: Transmition rate
            mu: E/I initial rate
            alpha: Movility function (quarantine object function)
            k: Saturation Kynetics Dynamics factor
            H_cap: Hospital capacity, either an int or a function(t)
            V_cap: VMI capacity, either an int or a function(t)



    """
    def __init__(self, config = None, inputdata=None):
        
        # Parameter that defines sigmoid's raising speed
        self.gw=20
        self.config = config

        # Later add the option to initialize it with a default cfg file. 
        if not config:
            print('Missing configuration file ')
            return None

        # ------------------------------- #
        #         Load parameters         #
        # ------------------------------- #
        if type(config) == dict:
            self.cfg = config    
        else:
            self.cfg = toml.load(config)
                
        # Model
        self.model = self.cfg['model']

        # Import fixed variables
        self.__dict__.update(self.cfg['parameters']['static'])
        
        self.tsim = self.t_end - self.t_init

        # Build functions
        for key in self.cfg['parameters']['dynamic']:
            if type(self.cfg['parameters']['dynamic'][key])==str:
                self.__dict__.update({key:cv19functions.build(self.cfg['parameters']['dynamic'][key])})
            else:
                # For inputing functions directly in the dictionary
                self.__dict__.update({key:self.cfg['parameters']['dynamic'][key]})
        
        
        # Ephemeris
        if 'ephemeris' in self.cfg:
            self.ephemeris = self.cfg['ephemeris']

        # Data:
        self.__dict__.update(self.cfg['data'])
        self.initdate = cv19timeutils.txt2Datetime(self.initdate) 

        if inputdata:
            self.inputdata = inputdata
            self.data = self.inputdata.data
            self.initdate = self.inputdata.initdate
            #if not type(self.initdate) == datetime.datetime:
            #    self.initdate = cv19timeutils.txt2Datetime(self.initdate)                

            #self.country = self.inputdata.country 
            self.state = self.inputdata.tstate 
            #self.county = self.inputdata.county 
            #self.healtservice = self.inputdata.healtservice 
            #self.loc_name = self.inputdata.loc_name             

        else:
            if self.datafile:
                # Falta construir funcion para cargar datos desde archivos, pero primero tengo que construir ese archivo 
                self.data = pd.DataFrame(self.cfg['data']['datafile'])
                self.inputdata = None
            elif self.importdata:                	            
                #self.inputdata = cv19data.ImportData(country=self.country,state=self.state,county=self.county,healthservice=self.healthservice,initdate=self.initdate,user = None,password = None)
                self.inputdata = cv19data.ImportData(tstate=self.state,initdate=self.initdate,user = None,password = None)
                self.inputdata.importdata()
                self.data = self.inputdata.data
            else: 
                print('No external data added')
                self.data = None
                self.inputdata = None

        # ------------------------------- #
        #       Initial conditions        #
        # ------------------------------- #
        # Revisar!!
        self.IC = self.cfg['initialconditions']
        for key in self.IC:
            if type(self.IC[key]) == str:
                # Crear error cuando no haya archivo de datos
                self.__dict__.update({key:self.data[self.IC[key]][0]})
            else:
                self.__dict__.update({key:self.IC[key]})



        self.setrelationalvalues()
        self.setequations()

        print('SEIR object created')

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
        self.S = self.N-self.E-self.I        
                    

    def setnewparams(self):
        self.setequations()
        self.setrelationalvalues()
        print('State parameters updated')


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
        self.dS=lambda t,S,E,I,R: self.S_f(t) - self.alpha(t)*self.beta(t)*S*(self.expinfection*E+I)/(self.N+self.k_I*I + self.k_R*R) + self.eta(t)*R
        
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
        self.dR=lambda t,I,R: self.R_f(t) + I/self.tI_R(t) - self.eta(t)*R

        # 6) Recovered per day
        self.dR_d=lambda t,I,R_d: self.R_f(t) + I/self.tI_R(t) - R_d

        # Recovered Accumulated
        #self.dR_ac =lambda t,I: self.pI_R/self.tI_R(t)*I 

        # 7) External Flux:
        self.dFlux = lambda t: self.S_f(t) + self.E_f(t) + self.I_f(t) + self.R_f(t) 


    """
    def calculateindicators(self):
        self.R_ef
        self.SHFR
        # Peak
        self.peak
        self.peakindex
        self.peak_t
        self.peak_date


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

    def integr(self,t0,T,h,E0init=False):
        self.integrate()
        return

    # Scipy
    def integrate(self,t0=0,T=None,h=0.01):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,mu*(self.I_ac)T] function doesn't
        #start. Or it's start if class object is initialze.

        if T is None:
            T = self.tsim

        if(not isinstance(self.S, np.ndarray)):

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
            return()
            

        
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

        initcond = np.array([S0,E0,E_d0,I0,I_d0,R0,R_d0,Flux0]) 
        
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

        self.anaytics()
        self.dfbuild(sol)

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

        if not isinstance(self.S, np.ndarray):

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

        self.anaytics()
        self.dfbuild(sol)
        return(sol)

    def anaytics(self):
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
        self.resume = pd.DataFrame({'peak':int(self.peak),'peak_t':self.peak_t,'peak_date':self.peak_date},index=[0])