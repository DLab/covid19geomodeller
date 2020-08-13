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

class simSEIRHVD:
    definputarray=np.array([
            [500,0.8,0.6,0,0,500,0],
            [500,0.6,0.5,0,0,500,0],
            [500,0.4,0.4,0,0,500,0]])            


    def __init__(self,beta = 0.15, mu =1,inputarray = definputarray,I=0,R=0,expinfection=1,SeroPrevFactor=1,population=100000,
    intgr = 0,I_ac =0, k =0):
        self.mu = mu
        self.beta = beta 
        self.sims = []
        self.inputarray=inputarray
        self.simulated = False
        self.I = I        
        self.R  = R    
        
        self.expinfection = expinfection 
        self.SeroPrevFactor=SeroPrevFactor
        self.population = population
        self.intgr = intgr

        # Accumulated Infected
        self.I_ac = I_ac

        # dayly Infected
        self.I_d = 0

        # Saturated Kinetics
        self.k = k

         
    
    def sim_run(self,tsim,max_mov,rem_mov,qp,iqt=0,fqt = 300,movfunct = 0):       
        case = SEIRHUDV(tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct,k=self.k)
        case.beta = self.beta       
        case.mu = self.mu
        case.I = self.I 
        case.R  = self.R 

        case.expinfection = self.expinfection
        case.SeroPrevFactor = self.SeroPrevFactor
        case.population = self.population 

        # Accumulated Infected
        case.I_ac = self.I_ac 
        
        case.I_d = self.I_d       

        #case.k = self.k
        
        case.setrelationalvalues()
        if self.intgr == 0:
            print('Fast Solver')
            case.integr_sci(0,tsim,0.1,False)        
        else:
            print('Robust Solver')            
            case.integr(0,tsim,0.1,False)
        out=[case,max_mov,rem_mov,qp,tsim]
        return(out)   
    
    def simulate(self,intgr=0):
        num_cores = multiprocessing.cpu_count()
        #params=Parallel(n_jobs=num_cores, verbose=50)(delayed(ref_test.refinepso_all)(Ir,tr,swarmsize=200,maxiter=50,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,1],Q_r=[0,1],obj_func='IN')for i in range(int(rep)))
        self.sims=Parallel(n_jobs=num_cores, verbose=50)(delayed(self.sim_run)(self.inputarray[i,0],self.inputarray[i,1],self.inputarray[i,2],self.inputarray[i,3],self.inputarray[i,4],self.inputarray[i,5],self.inputarray[i,6]) for i in range(self.inputarray.shape[0]))
        self.simulated = True
        return(self.sims)





class SEIR:  
    def __init__(self,tsim,alpha,beta,mu,k=0,I=100,I_ac=0,I_d=0,R=0,population=1000000,expinfection = 1, SeroPrevFactor=1):
        
        self.tsim = tsim
        self.alpha = alpha
        self.beta = beta
        self.mu=1.4
        self.k = k

        self.I = I                
        self.I_ac = I_ac        
        self.I_d = I_d

        self.R = R 
 
        self.SeroPrevFactor = SeroPrevFactor
        self.population = population        
        self.expinfection = expinfection         
        
        self.t=0

        # Expuestos
        self.E = self.mu*self.I               
        
        # Valores globales
        self.N =  self.SeroPrevFactor*self.population
        self.S = self.N- self.E-self.I - self.R               
        
        self.setparams()

        #self.setinitvalues()          
        #self.setscenario(tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct)                        
        
        
        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #        
        # dVariable/dt = sum(prob_i/in_time_i*in_State_i,i in in states) - sum(prob_i/out_time_i*out_State_i,i in in states) 
        
        # Susceptibles
        # dS/dt:
        self.dS=lambda t,S,E,I,R: -self.alpha(t)*self.beta*S*(self.expinfection*(E)+I)/(self.N+self.k*I)+self.eta*R
        # Exposed
        # dE/dt
        self.dE=lambda t,S,E,I: self.alpha(t)*self.beta*S*(self.expinfection*(E)+I)/(self.N+self.k*I)\
            -self.pasas/self.tasas*E

        # Infected
        # dI_as/dt
        self.dI=lambda t,E,I: self.pasas/self.tasas*E-self.pasR/self.tasR*I

        # Recovered
        # dR/dt
        self.dR=lambda t,I,R: self.pasR/self.tasR*I+-self.eta*R

        #Auxiliar functions:
        self.dI_ac=lambda t,E: self.pasas/self.tasas*E

        # Daily Infected
        self.dI_d = lambda t,E,I_d: self.pasas/self.tasas*E - I_d 
    


    def setparams(self):
        self.mu = 2.6

        self.beta = 0.19 # (*probabilidad de transmision por contacto con contagiados*)
        self.betaD = 0.0 #(*probabilidad de transmision por contacto con muertos*)

        self.pSas = 1#3 # Transicion de Susceptible a Expuesto Asintomatico
        self.tSas = 1.0

        self.pSsy = 0#7 # Transicion de Susceptible a Expuesto sintomatico
        self.tSsy = 1.0

        self.pasas = 1.0# Transicion de Expuesto asintomatico a Infectado asintomatico
        self.tasas = 5.0

        self.psymi = 0.8333#78 # Transicion de Expuesto Sintomatico a Infectado Mild
        self.tsymi = 5.0

        self.psycr = 0.05555#8 # Transicion de Expuesto Sintomatico a Infectado critico
        self.tsycr = 3.0

        self.psyse = 0.11111#14 # Transicion de Expuesto Sintomatico a Infectado Severo
        self.tsyse = 3.0

        self.pasR = 1.0   # Transicion de Infectado asintomatico a Recuperado
        self.tasR = 10.0 

        self.pmiR = 1.0  # Transicion de Infectado mild a Recuperado
        self.tmiR = 15.0

        self.psein = 1.0  # Transicion de Infectado serio a Hospitalizado (si no ha colapsado Htot)
        self.tsein = 3.0 

        self.pincrin = 0.03 # Transicion de Hospitalizado a Hospitalizado Critico (si no ha colapsado Htot)
        self.tincrin = 3.0

        self.pcrcrin = 1.0 # Transicion de Infectado critico  a Hopsitalizado Critico (si no ha colapsado Htot)
        self.tcrcrin = 3.0 

        self.pcrinV = 1.0 # Transicion de Hospitalizado critico a Ventilado (si no ha colapsado V)
        self.tcrinV = 0.01 

        self.pcrinD = 1.0 # Muerte de hospitalizado critico (Cuando V colapsa)
        self.tcrinD = 0.001 #

        self.pcrD = 1.0 # Muerte de Infectado critico (si ha colapsado Htot)
        self.tcrD = 3.0 #(*Hin+H_cr_in+Hout colapsa*)

        self.pseD = 1.0 # Muerte de Infectado severo (si ha colapsado Htot)
        self.tseD = 3.0

        self.pinout = 0.97 # Mejora de paciente severo hospitalizado, transita a Hout
        self.tinout = 6.0

        self.pVout = 0.5 # Mejora de ventilado hospitalizado, transita a Hout
        self.tVout = 15.0

        self.pVD = 0.5 # Muerte de ventilado
        self.tVD = 15.0

        self.poutR = 1.0 # Mejora del paciente hospitalizado, Hout a R
        self.toutR = 5.0

        self.pDB = 1.0 # Entierros
        self.tDB = 1.0 

        self.pcrV = 1.0 # Transicion de Hospitalizado critico a Ventilado (si no ha colapsado V)
        self.tcrV = 3.0 

        self.eta = 0.0 # tasa de perdida de inmunidad (1/periodo)


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

    def setscenario(self,tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct):
        self.tsim = tsim
        self.max_mov = max_mov
        self.rem_mov = rem_mov
        self.qp = qp
        self.iqt = iqt
        self.fqt = fqt
        if movfunct == 0:
            self.movfunct = 'once'
        elif movfunct == 1:
            self.movfunct = 'square'
        elif movfunct == 2:
            self.movfunct = 'sawtooth'
        else:
            self.movfunct = 'once'
        
        self.alpha = self.alphafunct(self.max_mov,self.rem_mov,self.qp,self.iqt,self.fqt,self.movfunct)
        return() 
        
    def alphafunct(self,max_mov,rem_mov,qp,iqt=0,fqt=300,movfunct = 'once'):
        """    
        # max_mov: Movilidad sin cuarentena
        # rem_mov: Movilidad con cuarentena
        # qp: Periodo cuarentena dinamica 
        #          - qp >0 periodo Qdinamica 
        #          - qp = 0 sin qdinamica
        # iqt: Initial quarantine time. Tiempo inicial antes de cuarentena dinamica
        #          - iqt>0 inicia con cuarentena total hasta iqt
        #          - iqt<0 sin cuarentena hasta iqt
        # fqt: Final quarantine time. Duracion tiempo cuarentena 
        # movfunct: Tipo de cuarentena dinamica desde iqt
        #          - once: una vez durante qp dias 
        #          - total: total desde iqt
        #          - sawtooth: diente de cierra
        #          - square: onda cuadrada
        """
        def alpha(t):             
            if 'square' in movfunct:
               def f(t): 
                   return signal.square(t)
               if t<abs(iqt):
                   if iqt>0:
                       return(rem_mov)
                   else:
                       return(max_mov)
               else:
                   if qp == 0:
                       return(max_mov)
                   elif t<fqt:
                       return((max_mov-rem_mov)/2*(f(np.pi / qp * t - np.pi))+(max_mov+rem_mov)/2)
                   else:
                       return(max_mov)   


            elif 'once' in movfunct:        
                if t<iqt:
                    return(max_mov)
                elif t>fqt:
                    return(max_mov)
                else:
                    return(rem_mov)


            elif 'sawtooth' in movfunct:
               def f(t): 
                   return signal.sawtooth(t)
               if t<abs(iqt):
                   if iqt>0:
                       return(rem_mov)
                   else:
                       return(max_mov)
               else:
                   if qp == 0:
                       return(max_mov)
                   elif t<fqt:
                       return((max_mov-rem_mov)/2*(f(np.pi / qp * t - np.pi))+(max_mov+rem_mov)/2)
                   else:
                       return(max_mov)   
                     
        return(alpha)


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
            
            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)

        else:
            return()

        
        def model_SEIR_graph(t,y,ydot):
            
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3])
            ydot[1]=self.dE(t,y[0],y[1],y[2])
            ydot[2]=self.dI(t,y[1],y[2])
            ydot[3]=self.dR(t,y[2],y[3])
            
            ydot[4]=self.dI_ac(t,y[1])
            ydot[5]=self.dI_d(t,y[1],y[5])

            
        initcond = np.array([S0,E0,I0,R0,I_ac0,I_d0])                                


        sol = odeint(model_SEIR_graph, self.t, initcond,method='admo')
        
        self.t=sol.values.t 
        
        self.S=sol.values.y[:,0]
        self.E=sol.values.y[:,1]        
        self.I=sol.values.y[:,2]
        self.R=sol.values.y[:,3]
        self.I_ac=sol.values.y[:,4]
        self.I_d=sol.values.y[:,5]                        
               
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

            self.t=np.arange(t0,T+h,h)

        else:
            return()
            

        
        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3])
            ydot[1]=self.dE(t,y[0],y[1],y[2])
            ydot[2]=self.dI(t,y[1],y[2])
            ydot[3]=self.dR(t,y[2],y[3])
            
            ydot[4]=self.dI_ac(t,y[1])
            ydot[5]=self.dI_d(t,y[1],y[5])
                         
                                          
            return(ydot)
        initcond = np.array([S0,E0,I0,R0,I_ac0,I_d0])  

        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA')
        
        self.t=sol.t 

         
        self.S=sol.y[0,:]
        self.E=sol.y[1,:]
        self.I=sol.y[2,:]
        self.R=sol.y[3,:]
        self.I_ac=sol.y[4,:]
        self.I_d=sol.y[5,:]

        return(sol)







