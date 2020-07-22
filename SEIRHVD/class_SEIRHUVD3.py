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


    def __init__(self,beta = 0.19, mu =2.6,inputarray = definputarray,B=221,D=26,V=758,I_act0=12642,R=0,Htot=None,Vtot=None,H_cr=80,H0=1720,tsat=30,Hmax=4000,Vmax=2000,expinfection=0,SeroPrevFactor=1,population=100000,
    intgr = 0,I_as_ac =0, I_mi_ac = 0, I_se_ac = 0, I_cr_ac = 0,H_crD = 0, VD=0,I_crD=0,I_seD=0,I_as_prop = 0.35, I_mi_prop = 0.63,I_cr_prop = 0.007,I_se_prop = 0.013, k =0):
        self.mu = mu
        self.beta = beta 
        self.sims = []
        self.inputarray=inputarray
        self.simulated = False
        self.B = B
        self.D = D
        self.V = V
        self.I_act0 = I_act0        
        self.R  = R    
        self.H_cr = H_cr # Hospitalizados criticos dia 0
        self.H0 = H0 # Hospitalizados totales dia 0
        
        self.Hmax = Hmax
        self.Vmax = Vmax
        self.expinfection = expinfection 
        self.SeroPrevFactor=SeroPrevFactor
        self.population = population
        self.Htot = Htot
        self.Vtot = Vtot
        self.intgr = intgr

        # Accumulated Infected
        self.I_as_ac = I_as_ac
        self.I_mi_ac = I_mi_ac
        self.I_se_ac = I_se_ac
        self.I_cr_ac = I_cr_ac

        # Deaths
        self.H_crD = H_crD
        self.VD = VD
        self.I_seD = I_seD
        self.I_crD = I_crD

        # Initial Infected proportion
        self.I_as_prop = I_as_prop
        self.I_mi_prop = I_mi_prop
        self.I_cr_prop = I_cr_prop
        self.I_se_prop = I_se_prop


        # dayly Infected
        self.I_as_d = 0
        self.I_mi_d = 0
        self.I_se_d = 0
        self.I_cr_d = 0

        # Saturated Kinetics
        self.k = k                        

         
    
    def sim_run(self,tsim,max_mov,rem_mov,qp,iqt=0,fqt = 300,movfunct = 0):       
        case = SEIRHUDV(tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct,k=self.k)
        case.beta = self.beta       
        case.mu = self.mu
        case.B = self.B 
        case.D = self.D 
        case.V = self.V 
        case.I_act0 = self.I_act0 
        case.R  = self.R 
        case.Htot = self.Htot
        case.Vtot = self.Vtot         
        case.H_cr = self.H_cr  # Hospitalizados criticos dia 0
        case.H0 = self.H0  # Hospitalizados totales dia 0
        
        case.Hmax = self.Hmax
        case.Vmax = self.Vmax
        case.expinfection = self.expinfection
        case.SeroPrevFactor = self.SeroPrevFactor
        case.population = self.population 

        # Accumulated Infected
        case.I_as_ac = self.I_as_ac 
        case.I_mi_ac = self.I_mi_ac 
        case.I_se_ac = self.I_se_ac 
        case.I_cr_ac = self.I_cr_ac 

        # Deaths
        case.H_crD = self.H_crD
        case.VD = self.VD
        case.I_seD = self.I_seD
        case.I_crD = self.I_crD

        case.I_as_d = self.I_as_d
        case.I_mi_d = self.I_mi_d
        case.I_se_d = self.I_se_d
        case.I_cr_d = self.I_cr_d

        case.I_as_prop=self.I_as_prop
        case.I_mi_prop=self.I_mi_prop
        case.I_cr_prop=self.I_cr_prop
        case.I_se_prop=self.I_se_prop
        

        case.setrelationalvalues()
        if self.intgr == 0:
            print('Scipy Solver')
            case.integr_sci(0,tsim,0.1,False)        
        else:
            print('Scikits-odes solver (robust)')                  
            case.integr(0,tsim,0.1,False)
        out=[case,max_mov,rem_mov,qp,tsim]
        return(out)   


        # Agregar un codigo numerico para que solo haya 1 vector de simulacion
        #['once','once','once','once','once','once','square','square','once','once','once','once','once','once','square','square']
        
    def setparameters(self):
        return
    def simulate(self,intgr=0):
        num_cores = multiprocessing.cpu_count()
        #params=Parallel(n_jobs=num_cores, verbose=50)(delayed(ref_test.refinepso_all)(Ir,tr,swarmsize=200,maxiter=50,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,1],Q_r=[0,1],obj_func='IN')for i in range(int(rep)))
        self.sims=Parallel(n_jobs=num_cores, verbose=50)(delayed(self.sim_run)(self.inputarray[i,0],self.inputarray[i,1],self.inputarray[i,2],self.inputarray[i,3],self.inputarray[i,4],self.inputarray[i,5],self.inputarray[i,6]) for i in range(self.inputarray.shape[0]))
        self.simulated = True
        return(self.sims)

    def getscenarios(self):
        return()
    def addscenario(self,inputarray):
        return()





class SEIRHUDV :  
    def __init__(self,tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct,k=0):
        #print(max_mov)
        self.setparams()
        self.setinitvalues()  
        #print(max_mov)
        self.setscenario(tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct)
        self.k=k        
        #global mobility reduction parameters
        #self.alpha=alpha
        # total Hostipatl beds

        
        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #
        
        # dVariable/dt = sum(prob_i/in_time_i*in_State_i,i in in states) - sum(prob_i/out_time_i*out_State_i,i in in states) 
        
        # Susceptibles
        # dS/dt:
        self.dS=lambda t,S,E_as,E_sy,I_as,I_mi,D,R: -self.alpha(t)*self.beta*S*(self.expinfection*(E_as+E_sy)+I_as+I_mi)/(self.N+self.k*(I_as+I_mi))-self.betaD*D+self.eta*R
        
        # Exposed
        # dE_as/dt
        self.dE_as=lambda t,S,E_as,E_sy,I_as,I_mi: self.pSas/self.tSas*self.alpha(t)*self.beta*S*(self.expinfection*(E_as+E_sy)+I_as+I_mi)/(self.N+self.k*(I_as+I_mi))\
            -self.pasas/self.tasas*E_as
        # dE_sy/dt
        self.dE_sy=lambda t,S,E_as,E_sy,I_as,I_mi: self.pSsy/self.tSsy*self.alpha(t)*self.beta*S*(self.expinfection*(E_as+E_sy)+I_as+I_mi)/(self.N+self.k*(I_as+I_mi))\
            -self.psymi/self.tsymi*E_sy-self.psyse/self.tsyse*E_sy-self.psycr/self.tsycr*E_sy
        
        # Infected
        # dI_as/dt
        self.dI_as=lambda t,E_as,I_as: self.pasas/self.tasas*E_as-self.pasR/self.tasR*I_as
        # dI_mi/dt
        self.dI_mi=lambda t,E_sy,I_mi: self.psymi/self.tsymi*E_sy-self.pmiR/self.tmiR*I_mi
        # dI_se/dt: Esy -  
        self.dI_se=lambda t,E_sy,I_se,H_in,H_cr,H_out: self.psyse/self.tsyse*E_sy-self.psein/self.tsein*I_se*(self.h_sat(H_in,H_cr,H_out,t))\
            -self.pseD/self.tseD*I_se*(1-self.h_sat(H_in,H_cr,H_out,t))            
        # dI_cr/dt
        self.dI_cr=lambda t,E_sy,I_cr,H_in,H_cr,H_out,V: self.psycr/self.tsycr*E_sy\
            -self.pcrcrin/self.tcrcrin*I_cr*(self.h_sat(H_in,H_cr,H_out,t)*(1-self.v_sat(V,t)))\
            -self.pcrV/self.tcrV*I_cr*(self.v_sat(V,t))-self.pcrD/self.tcrD*I_cr*(1-self.h_sat(H_in,H_cr,H_out,t))*(1-self.v_sat(V,t))

        #self.dI_cr=lambda t,E_sy,I_cr,H_in,H_cr,H_out: self.psycr/self.tsycr*E_sy-self.pcrcrin/self.tcrcrin*I_cr*(self.h_sat(H_in,H_cr,H_out,t))\
        #    -self.pcrD/self.tcrD*I_cr*(1-self.h_sat(H_in,H_cr,H_out,t))
        
        # Hospitalized
        # dH_in/dt: Hospitalized arriving to the hospital
        self.dH_in=lambda t,I_se,H_in,H_cr,H_out: self.psein/self.tsein*I_se*(self.h_sat(H_in,H_cr,H_out,t))\
            -self.pinout/self.tinout*H_in-self.pincrin/self.tincrin*H_in
        # dH_cr/dt: Hospitalized in critical conditions            
        self.dH_cr=lambda t,I_cr,H_in,H_cr,H_out,V: self.pcrcrin/self.tcrcrin*I_cr*(self.h_sat(H_in,H_cr,H_out,t)*(1-self.v_sat(V,t)))\
            -self.pcrinV/self.tcrinV*H_cr*self.v_sat(V,t)-self.pcrinD/self.tcrinD*H_cr*(1-self.v_sat(V,t))+self.pincrin/self.tincrin*H_in
        # dH_out/dt: Hospitalized getting better
        self.dH_out=lambda t,H_in,H_out,V: self.pVout/self.tVout*V+self.pinout/self.tinout*H_in-self.poutR/self.toutR*H_out
        # dV/dt: Ventilator necessity rates
        self.dV=lambda t,I_cr,H_cr,V: self.pcrinV/self.tcrinV*H_cr*self.v_sat(V,t)-self.pVout/self.tVout*V-self.pVD/self.tVD*V\
        +self.pcrV/self.tcrV*I_cr*(self.v_sat(V,t))
        
        # Deaths
        # dD/dt: Death Rate
        self.dD=lambda t,I_se,I_cr,H_in,H_cr,H_out,V,D: self.pVD/self.tVD*V+self.pcrD/self.tcrD*I_cr*(1-self.h_sat(H_in,H_cr,H_out,t))*(1-self.v_sat(V,t))+\
            self.pseD/self.tseD*I_se*(1-self.h_sat(H_in,H_cr,H_out,t))+self.pcrinD/self.tcrinD*H_cr*(1-self.v_sat(V,t))-\
            self.pDB/self.tDB*D
        # dB/dt: Bury rate
        self.dB=lambda t,D: self.pDB/self.tDB*D
        
        # Recovered
        # dR/dt
        self.dR=lambda t,I_as,I_mi,H_out,R: self.pasR/self.tasR*I_as+self.pmiR/self.tmiR*I_mi+self.poutR/self.toutR*H_out-self.eta*R

        #Auxiliar functions:
        self.dI=lambda t,E_as,E_sy: self.pasas/self.tasas*E_as+self.psymi/self.tsymi*E_sy+self.psyse/self.tsyse*E_sy+self.psycr/self.tsycr*E_sy
        self.dCV=lambda t,I_cr,H_in,H_cr,H_out,V,CV:self.pcrD/self.tcrD*I_cr*(1-self.h_sat(H_in,H_cr,H_out,t)*(1-self.v_sat(V,t)))+self.pcrinD/self.tcrinD*H_cr*(1-self.v_sat(V,t))-CV
        self.dCH=lambda t,I_se,H_in,H_cr,H_out,CH:self.pseD/self.tseD*I_se*(1-self.h_sat(H_in,H_cr,H_out,t))-CH
        self.dACV=lambda t,CV: CV
        self.dACH=lambda t,CH: CH

        # Accummulated Infected
        self.dI_as_ac = lambda t,E_as: self.pasas/self.tasas*E_as
        self.dI_mi_ac = lambda t,E_sy: self.psymi/self.tsymi*E_sy
        self.dI_se_ac = lambda t,E_sy: self.psyse/self.tsyse*E_sy
        self.dI_cr_ac = lambda t,E_sy: self.psycr/self.tsycr*E_sy

        # Accumulated Deads         
        self.dH_crD=lambda t,H_cr,V: self.pcrinD/self.tcrinD*H_cr*(1-self.v_sat(V,t)) 
        self.dVD=lambda t,V: self.pVD/self.tVD*V        
        self.dI_seD=lambda t,I_se,H_in,H_cr,H_out: self.pseD/self.tseD*I_se*(1-self.h_sat(H_in,H_cr,H_out,t))
        self.dI_crD=lambda t,I_cr,H_in,H_cr,H_out,V: self.pcrD/self.tcrD*I_cr*(1-self.h_sat(H_in,H_cr,H_out,t))*(1-self.v_sat(V,t))



        # Daily Infected
        self.dI_as_d = lambda t,E_as,I_as_d: self.pasas/self.tasas*E_as - I_as_d 
        self.dI_mi_d = lambda t,E_sy,I_mi_d: self.psymi/self.tsymi*E_sy - I_mi_d 
        self.dI_se_d = lambda t,E_sy,I_se_d: self.psyse/self.tsyse*E_sy - I_se_d
        self.dI_cr_d = lambda t,E_sy,I_cr_d: self.psycr/self.tsycr*E_sy - I_cr_d

        # Daily Deads         
        self.dH_crD_d=lambda t,H_cr,V,H_crD_d: self.pcrinD/self.tcrinD*H_cr*(1-self.v_sat(V,t)) - H_crD_d
        self.dVD_d=lambda t,V,VD_d: self.pVD/self.tVD*V - VD_d
        self.dI_seD_d=lambda t,I_se,H_in,H_cr,H_out,I_seD_d: self.pseD/self.tseD*I_se*(1-self.h_sat(H_in,H_cr,H_out,t)) - I_seD_d
        self.dI_crD_d=lambda t,I_cr,H_in,H_cr,H_out,V,I_crD_d: self.pcrD/self.tcrD*I_cr*(1-self.h_sat(H_in,H_cr,H_out,t))*(1-self.v_sat(V,t)) - I_crD_d

        


    # UCI and UTI beds saturation function
    def h_sat(self,H_in,H_cr,H_out,t):
        return(expit(-self.gw*(H_in+H_cr+H_out-self.Htot(t))))
    # Ventilators Saturation Function    
    def v_sat(self,V,t):
        return(expit(-self.gw*(V-self.Vtot(t))))

    def setparams(self):
        self.mu = 2.6

        self.beta = 0.19 # (*probabilidad de transmision por contacto con contagiados*)
        self.betaD = 0.0 #(*probabilidad de transmision por contacto con muertos*)

        self.pSas = 0.4#3 # Transicion de Susceptible a Expuesto Asintomatico
        self.tSas = 1.0

        self.pSsy = 0.6#7 # Transicion de Susceptible a Expuesto sintomatico
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
        # 15 de Mayo
        self.I_act0 = 12642
        self.res=1 
        self.muS=self.mu        
        self.H_incr = 34
        self.H_incr2 = 0
        self.H_incr3 = 0
        self.V_incr = 17
        self.V_incr2 = 0
        self.V_incr3 = 0
        self.Vc0 = 1029
        self.Hc0 = 1980
        self.H0=1720 #1980#1903.0
        self.H_cr=80.0
        self.gw=10
        self.D=26.0
        self.B=221.0
        self.R=0.0
        self.V=758.0#846.0
        self.mu=1.4
        self.t=400.0
        self.CV=0
        self.CH=0
        self.ACV=0
        self.ACH=0
        self.SeroPrevFactor = 1
        self.population = 8125072
        
        self.Hmax = 3000
        self.Vmax = 1500
        self.expinfection = 0
        
        # Accumulated Infected
        self.I_as_ac = 0
        self.I_mi_ac = 0
        self.I_se_ac = 0
        self.I_cr_ac = 0

        # Deaths
        self.H_crD = 0
        self.VD = 0
        self.I_seD = 0
        self.I_crD = 0

        # Daily Infected
        self.I_as_d = 0
        self.I_mi_d = 0
        self.I_se_d = 0
        self.I_cr_d = 0

        # Daily Deaths
        self.H_crD_d = 0
        self.VD_d = 0
        self.I_seD_d = 0
        self.I_crD_d = 0

        # Initial Infected proportion
        self.I_as_prop = 0.35
        self.I_mi_prop = 0.63
        self.I_cr_prop = 0.007
        self.I_se_prop = 0.013

        self.setrelationalvalues()

    def setrelationalvalues(self):       
        self.I_as= self.I_as_prop*self.I_act0 
        self.I_mi= self.I_mi_prop*self.I_act0 
        self.I_cr= self.I_cr_prop*self.I_act0 
        self.I_se = self.I_se_prop*self.I_act0              
        # Expuestos

        self.E_as=0.3*self.mu*self.I_act0
        self.E_sy=0.7*self.mu*self.I_act0 

        # Hospitalizados
        #self.V+=(self.I_act0-(self.I_as+self.I_mi+self.I_cr+self.I_se))*0.05
        self.H_in=self.H0*0.42-self.H_cr/2 #+ (self.I_act0-(self.I_as+self.I_mi+self.I_cr+self.I_se))*0.1
        self.H_out=self.H0*0.58-self.H_cr/2 
       
        # Valores globales
        self.SeroPrevPop =  self.SeroPrevFactor*self.population
        self.S=self.SeroPrevPop-self.H0-self.V-self.D-(self.E_as+self.E_sy)-(self.I_as+self.I_cr+self.I_se+self.I_mi)
        self.N=(self.S+self.E_as+self.E_sy+self.I_as+self.I_mi+self.I_se+self.I_cr+self.H_in+self.H_cr+self.H_out+self.V+self.D+self.R)
        self.I=self.I_cr+self.I_as+self.I_se+self.I_mi        

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



    def get_conditions(self):
        # This function will return a text explaining the different simulation scenarios currently placed as inputs            
        return

 
    def integr(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.
        print('Import odeint')
        from scikits.odes.odeint import odeint


        if(not isinstance(self.S, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                E_as0=self.mu*(self.I_as)
                E_sy0=self.mu*(self.I_mi+self.I_se+self.I_cr)
            else:
                E_as0=self.E_as
                E_sy0=self.E_sy
            S0=self.S
            I_as0=self.I_as
            I_mi0=self.I_mi
            I_se0=self.I_se
            I_cr0=self.I_cr
            H_in0=self.H_in
            H_cr0=self.H_cr
            H_out0=self.H_out
            V0=self.V
            D0=self.D
            B0=self.B
            R0=self.R
            I0=self.I
            CV0=self.CV
            CH0=self.CH
            ACV0=self.ACV
            ACH0=self.ACH

            I_as_ac0=self.I_as_ac
            I_mi_ac0=self.I_mi_ac
            I_se_ac0=self.I_se_ac
            I_cr_ac0=self.I_cr_ac

            H_crD0 = self.H_crD
            VD0 = self.VD
            I_crD0 = self.I_crD
            I_seD0 = self.I_seD

            I_as_d0=self.I_as_d
            I_mi_d0=self.I_mi_d
            I_se_d0=self.I_se_d
            I_cr_d0=self.I_cr_d

            H_crD_d0 = self.H_crD_d
            VD_d0 = self.VD_d
            I_crD_d0 = self.I_crD_d
            I_seD_d0 = self.I_seD_d

            self.t=np.arange(t0,T+h,h)
            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            S0=self.S[idx]
            E_as0=self.E_as
            E_sy0=self.E_sy
            I_as0=self.I_as[idx]
            I_mi0=self.I_mi[idx]
            I_se0=self.I_se[idx]
            I_cr0=self.I_cr[idx]
            H_in0=self.H_in[idx]
            H_cr0=self.H_cr[idx]
            H_out0=self.H_out[idx]
            V0=self.V[idx]
            D0=self.D[idx]
            B0=self.B[idx]
            R0=self.R[idx]
            I0=self.I[idx]
            CV0=self.CV[idx]
            CH0=self.CH[idx]
            ACV0=self.ACV[idx]
            ACH0=self.ACH[idx]

            I_as_ac0=self.I_as_ac[idx]
            I_mi_ac0=self.I_mi_ac[idx]
            I_se_ac0=self.I_se_ac[idx]
            I_cr_ac0=self.I_cr_ac[idx]

            H_crD0 = self.H_crD[idx]
            VD0 = self.VD[idx]
            I_seD0 = self.I_seD[idx]
            I_crD0 = self.I_crD[idx]

            I_as_d0=self.I_as_d[idx]
            I_mi_d0=self.I_mi_d[idx]
            I_se_d0=self.I_se_d[idx]
            I_cr_d0=self.I_cr_d[idx]

            H_crD_d0 = self.H_crD_d[idx]
            VD_d0 = self.VD_d[idx]
            I_seD_d0 = self.I_seD_d[idx]
            I_crD_d0 = self.I_crD_d[idx]                                               

            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)

        else:
            return()

        
        def model_SEIR_graph(t,y,ydot):
            
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3],y[4],y[11],y[13])
            ydot[1]=self.dE_as(t,y[0],y[1],y[2],y[3],y[4])
            ydot[2]=self.dE_sy(t,y[0],y[1],y[2],y[3],y[4])
            ydot[3]=self.dI_as(t,y[1],y[3])
            ydot[4]=self.dI_mi(t,y[2],y[4])
            ydot[5]=self.dI_se(t,y[2],y[5],y[7],y[8],y[9])
            ydot[6]=self.dI_cr(t,y[2],y[6],y[7],y[8],y[9],y[10]) 
            ydot[7]=self.dH_in(t,y[5],y[7],y[8],y[9])
            ydot[8]=self.dH_cr(t,y[6],y[7],y[8],y[9],y[10])
            ydot[9]=self.dH_out(t,y[7],y[9],y[10])
            ydot[10]=self.dV(t,y[6],y[8],y[10])
            ydot[11]=self.dD(t,y[5],y[6],y[7],y[8],y[9],y[10],y[11])
            ydot[12]=self.dB(t,y[11])
            ydot[13]=self.dR(t,y[3],y[4],y[9],y[13])
            ydot[14]=self.dI(t,y[1],y[2])
            ydot[15]=self.dCV(t,y[6],y[7],y[8],y[9],y[10],y[15])
            ydot[16]=self.dCH(t,y[5],y[7],y[8],y[9],y[16])
            ydot[17]=self.dACH(t,y[15])
            ydot[18]=self.dACV(t,y[16])

            ydot[19]=self.dI_as_ac(t,y[1])
            ydot[20]=self.dI_mi_ac(t,y[2]) 
            ydot[21]=self.dI_se_ac(t,y[2]) 
            ydot[22]=self.dI_cr_ac(t,y[2])

            ydot[23]=self.dH_crD(t,y[8],y[10])
            ydot[24]=self.dVD(t,y[10])
            ydot[25]=self.dI_seD(t,y[5],y[7],y[8],y[9])
            ydot[26]=self.dI_crD(t,y[6],y[7],y[8],y[9],y[10])

            ydot[27]=self.dI_as_d(t,y[1],y[27])
            ydot[28]=self.dI_mi_d(t,y[2],y[28]) 
            ydot[29]=self.dI_se_d(t,y[2],y[29]) 
            ydot[30]=self.dI_cr_d(t,y[2],y[30])

            ydot[31]=self.dH_crD_d(t,y[8],y[10],y[31])
            ydot[32]=self.dVD_d(t,y[10],y[32])
            ydot[33]=self.dI_seD_d(t,y[5],y[7],y[8],y[9],y[33])
            ydot[34]=self.dI_crD_d(t,y[6],y[7],y[8],y[9],y[10],y[34])

            
        initcond = np.array([S0,E_as0,E_sy0,I_as0,I_mi0,I_se0,I_cr0,
                                H_in0,H_cr0,H_out0,V0,D0,B0,R0,I0,CV0,CH0,ACV0,ACH0,
                                I_as_ac0,I_mi_ac0,I_se_ac0,I_cr_ac0,H_crD0,VD0,I_seD0,I_crD0,I_as_d0,I_mi_d0,I_se_d0,I_cr_d0,H_crD_d0,VD_d0,I_seD_d0,I_crD_d0]) 


        sol = odeint(model_SEIR_graph, self.t, initcond,method='admo')
        
        self.t=sol.values.t 
        
        self.S=sol.values.y[:,0]
        self.E_as=sol.values.y[:,1]
        self.E_sy=sol.values.y[:,2]
        self.I_as=sol.values.y[:,3]
        self.I_mi=sol.values.y[:,4]
        self.I_se=sol.values.y[:,5]
        self.I_cr=sol.values.y[:,6]                
        self.H_in=sol.values.y[:,7]
        self.H_cr=sol.values.y[:,8]
        self.H_out=sol.values.y[:,9]
        self.V=sol.values.y[:,10]
        self.D=sol.values.y[:,11]
        self.B=sol.values.y[:,12]
        self.R=sol.values.y[:,13]
        self.I=sol.values.y[:,14]
        self.CV=sol.values.y[:,15]
        self.CH=sol.values.y[:,16]
        self.ACV=sol.values.y[:,17]
        self.ACH=sol.values.y[:,18]

        self.I_as_ac=sol.values.y[:,19]
        self.I_mi_ac=sol.values.y[:,20]
        self.I_se_ac=sol.values.y[:,21]
        self.I_cr_ac=sol.values.y[:,22]

        self.H_crD = sol.values.y[:,23]
        self.VD = sol.values.y[:,24]
        self.I_seD = sol.values.y[:,25]
        self.I_crD = sol.values.y[:,26]

        self.I_as_d=sol.values.y[:,27]
        self.I_mi_d=sol.values.y[:,28]
        self.I_se_d=sol.values.y[:,29]
        self.I_cr_d=sol.values.y[:,30]

        self.H_crD_d = sol.values.y[:,31]
        self.VD_d = sol.values.y[:,32]
        self.I_seD_d = sol.values.y[:,33]
        self.I_crD_d = sol.values.y[:,34]                
               
        return(sol)

    def integr_sci(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.

        if(not isinstance(self.S, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                E_as0=self.mu*(self.I_as)
                E_sy0=self.mu*(self.I_mi+self.I_se+self.I_cr)
            else:
                E_as0=self.E_as
                E_sy0=self.E_sy
            S0=self.S
            I_as0=self.I_as
            I_mi0=self.I_mi
            I_se0=self.I_se
            I_cr0=self.I_cr
            H_in0=self.H_in
            H_cr0=self.H_cr
            H_out0=self.H_out
            V0=self.V
            D0=self.D
            B0=self.B
            R0=self.R
            I0=self.I
            CV0=self.CV
            CH0=self.CH
            ACV0=self.ACV
            ACH0=self.ACH

            I_as_ac0=self.I_as_ac
            I_mi_ac0=self.I_mi_ac
            I_se_ac0=self.I_se_ac
            I_cr_ac0=self.I_cr_ac

            H_crD0 = self.H_crD
            VD0 = self.VD
            I_crD0 = self.I_crD
            I_seD0 = self.I_seD

            I_as_d0=self.I_as_d
            I_mi_d0=self.I_mi_d
            I_se_d0=self.I_se_d
            I_cr_d0=self.I_cr_d

            H_crD_d0 = self.H_crD_d
            VD_d0 = self.VD_d
            I_crD_d0 = self.I_crD_d
            I_seD_d0 = self.I_seD_d                          

            self.t=np.arange(t0,T+h,h)
            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            S0=self.S[idx]
            E_as0=self.E_as
            E_sy0=self.E_sy
            I_as0=self.I_as[idx]
            I_mi0=self.I_mi[idx]
            I_se0=self.I_se[idx]
            I_cr0=self.I_cr[idx]
            H_in0=self.H_in[idx]
            H_cr0=self.H_cr[idx]
            H_out0=self.H_out[idx]
            V0=self.V[idx]
            D0=self.D[idx]
            B0=self.B[idx]
            R0=self.R[idx]
            I0=self.I[idx]
            CV0=self.CV[idx]
            CH0=self.CH[idx]                
            ACV0=self.ACV[idx]
            ACH0=self.ACH[idx]

            I_as_ac0=self.I_as_ac[idx]
            I_mi_ac0=self.I_mi_ac[idx]
            I_se_ac0=self.I_se_ac[idx]
            I_cr_ac0=self.I_cr_ac[idx]

            H_crD0 = self.H_crD[idx]
            VD0 = self.VD[idx]
            I_seD0 = self.I_seD[idx]
            I_crD0 = self.I_crD[idx]

            I_as_d0=self.I_as_d[idx]
            I_mi_d0=self.I_mi_d[idx]
            I_se_d0=self.I_se_d[idx]
            I_cr_d0=self.I_cr_d[idx]

            H_crD_d0 = self.H_crD_d[idx]
            VD_d0 = self.VD_d[idx]
            I_seD_d0 = self.I_seD_d[idx]
            I_crD_d0 = self.I_crD_d[idx]    

            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)


        else:
            return()
            

        
        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3],y[4],y[11],y[13])
            ydot[1]=self.dE_as(t,y[0],y[1],y[2],y[3],y[4])
            ydot[2]=self.dE_sy(t,y[0],y[1],y[2],y[3],y[4])
            ydot[3]=self.dI_as(t,y[1],y[3])
            ydot[4]=self.dI_mi(t,y[2],y[4])
            ydot[5]=self.dI_se(t,y[2],y[5],y[7],y[8],y[9])
            ydot[6]=self.dI_cr(t,y[2],y[6],y[7],y[8],y[9],y[10]) 
            ydot[7]=self.dH_in(t,y[5],y[7],y[8],y[9])
            ydot[8]=self.dH_cr(t,y[6],y[7],y[8],y[9],y[10])
            ydot[9]=self.dH_out(t,y[7],y[9],y[10])
            ydot[10]=self.dV(t,y[6],y[8],y[10])
            ydot[11]=self.dD(t,y[5],y[6],y[7],y[8],y[9],y[10],y[11])
            ydot[12]=self.dB(t,y[11])
            ydot[13]=self.dR(t,y[3],y[4],y[9],y[13])
            ydot[14]=self.dI(t,y[1],y[2])
            ydot[15]=self.dCV(t,y[6],y[7],y[8],y[9],y[10],y[15])
            ydot[16]=self.dCH(t,y[5],y[7],y[8],y[9],y[16])
            ydot[17]=self.dACH(t,y[15])
            ydot[18]=self.dACV(t,y[16])

            ydot[19]=self.dI_as_ac(t,y[1])
            ydot[20]=self.dI_mi_ac(t,y[2]) 
            ydot[21]=self.dI_se_ac(t,y[2]) 
            ydot[22]=self.dI_cr_ac(t,y[2])

            ydot[23]=self.dH_crD(t,y[8],y[10])
            ydot[24]=self.dVD(t,y[10])
            ydot[25]=self.dI_seD(t,y[5],y[7],y[8],y[9])
            ydot[26]=self.dI_crD(t,y[6],y[7],y[8],y[9],y[10])

            ydot[27]=self.dI_as_d(t,y[1],y[27])
            ydot[28]=self.dI_mi_d(t,y[2],y[28]) 
            ydot[29]=self.dI_se_d(t,y[2],y[29]) 
            ydot[30]=self.dI_cr_d(t,y[2],y[30])

            ydot[31]=self.dH_crD_d(t,y[8],y[10],y[31])
            ydot[32]=self.dVD_d(t,y[10],y[32])
            ydot[33]=self.dI_seD_d(t,y[5],y[7],y[8],y[9],y[33])
            ydot[34]=self.dI_crD_d(t,y[6],y[7],y[8],y[9],y[10],y[34])                             
                                          
            return(ydot)
        initcond = np.array([S0,E_as0,E_sy0,I_as0,I_mi0,I_se0,I_cr0,
                                H_in0,H_cr0,H_out0,V0,D0,B0,R0,I0,CV0,CH0,ACV0,ACH0,
                                I_as_ac0,I_mi_ac0,I_se_ac0,I_cr_ac0,H_crD0,VD0,I_crD0,I_seD0,I_as_d0,I_mi_d0,I_se_d0,I_cr_d0,H_crD_d0,VD_d0,I_seD_d0,I_crD_d0])

        
        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA')
        
        self.t=sol.t 
        
        self.S=sol.y[0,:]
        self.E_as=sol.y[1,:]
        self.E_sy=sol.y[2,:]
        self.I_as=sol.y[3,:]
        self.I_mi=sol.y[4,:]
        self.I_se=sol.y[5,:]
        self.I_cr=sol.y[6,:]
        self.H_in=sol.y[7,:]
        self.H_cr=sol.y[8,:]
        self.H_out=sol.y[9,:]
        self.V=sol.y[10,:]
        self.D=sol.y[11,:]
        self.B=sol.y[12,:]
        self.R=sol.y[13,:]
        self.I=sol.y[14,:]
        self.CV=sol.y[15,:]
        self.CH=sol.y[16,:]
        self.ACV=sol.y[17,:]
        self.ACH=sol.y[18,:]
        self.I_as_ac=sol.y[19,:]
        self.I_mi_ac=sol.y[20,:]
        self.I_se_ac=sol.y[21,:]
        self.I_cr_ac=sol.y[22,:]
        self.H_crD = sol.y[23,:]
        self.VD = sol.y[24,:]
        self.I_seD = sol.y[25,:]
        self.I_crD = sol.y[26,:] 
        self.I_as_d=sol.y[27,:]
        self.I_mi_d=sol.y[28,:]
        self.I_se_d=sol.y[29,:]
        self.I_cr_d=sol.y[30,:]   
        self.H_crD_d = sol.y[31,:]
        self.VD_d = sol.y[32,:]
        self.I_seD_d = sol.y[33,:]
        self.I_crD_d = sol.y[34,:]         

        return(sol)






















    # def integr_SciML(self,t0,T,h,E0init=False):
    #     #integrator function that star form t0 and finish with T with h as
    #     #timestep. If there aren't inital values in [t0,T] function doesn't
    #     #start. Or it's start if class object is initialze.


    #     if(not isinstance(self.S, np.ndarray)):
    #         #pass if object is initalized
    #         if(E0init):
    #             E_as0=self.mu*(self.I_as)
    #             E_sy0=self.mu*(self.I_mi+self.I_se+self.I_cr)
    #         else:
    #             E_as0=self.E_as
    #             E_sy0=self.E_sy
    #         S0=self.S
    #         I_as0=self.I_as
    #         I_mi0=self.I_mi
    #         I_se0=self.I_se
    #         I_cr0=self.I_cr
    #         H_in0=self.H_in
    #         H_cr0=self.H_cr
    #         H_out0=self.H_out
    #         V0=self.V
    #         D0=self.D
    #         B0=self.B
    #         R0=self.R
    #         print(self.H_out)
    #         self.t=np.arange(t0,T+h,h)
            
    #     elif((min(self.t)<=t0) & (t0<=max(self.t))):
    #         #Condition over exiting time in already initialized object

    #         #Search fot initial time
    #         idx=np.searchsorted(self.t,t0)

    #         #set initial condition

    #         S0=self.S[idx]
    #         E_as0=self.E_as
    #         E_sy0=self.E_sy
    #         I_as0=self.I_as[idx]
    #         I_mi0=self.I_mi[idx]
    #         I_se0=self.I_se[idx]
    #         I_cr0=self.I_cr[idx]
    #         H_in0=self.H_in[idx]

    #         #set time grid
    #         self.t=np.arange(self.t[idx],T+h,h)


    #     else:
    #         return()
            
    #     # dim=self.S.shape[0]
        
    #     def model_SEIR_graph(t,y):
    #         ydot=np.zeros(len(y))
    #         ydot[0]=self.dS(t,y[0],y[1],y[2],y[3],y[4],y[11],y[13])
    #         ydot[1]=self.dE_as(t,y[0],y[1],y[2],y[3],y[4])
    #         ydot[2]=self.dE_sy(t,y[0],y[1],y[2],y[3],y[4])
    #         ydot[3]=self.dI_as(t,y[1],y[3])
    #         ydot[4]=self.dI_mi(t,y[2],y[4])
    #         ydot[5]=self.dI_se(t,y[2],y[5],y[7],y[8],y[9])
    #         ydot[6]=self.dI_cr(t,y[2],y[6],y[7],y[8],y[9]) 
    #         ydot[7]=self.dH_in(t,y[5],y[7],y[8],y[9])
    #         ydot[8]=self.dH_cr(t,y[6],y[7],y[8],y[9],y[10])
    #         ydot[9]=self.dH_out(t,y[7],y[9],y[10])
    #         ydot[10]=self.dV(t,y[8],y[10])
    #         ydot[11]=self.dD(t,y[5],y[6],y[8],y[10],y[11])
    #         ydot[12]=self.dB(t,y[11])
    #         ydot[13]=self.dR(t,y[3],y[4],y[9],y[13])
    #         return(ydot)
    #     initcond = np.array([S0,E_as0,E_sy0,I_as0,I_mi0,I_se0,I_cr0,
    #                          H_in0,H_cr0,H_out0,V0,D0,B0,R0])

    #     # initcond = initcond.reshape(4*dim)
    #     numba_f = numba.jit(model_SEIR_graph)
    #     tspan=(float(t0),float(T))
    #     print(tspan)
    #     prob = de.ODEProblem(numba_f, initcond,(t0,T))
    #     sol = de.solve(prob,de.TRBDF2())
    #     # self.t=sol.t 

    #     # self.S=sol.y[0,:]
    #     # self.E_as=sol.y[1,:]
    #     # self.E_sy=sol.y[2,:]
    #     # self.I_as=sol.y[3,:]
    #     # self.I_mi=sol.y[4,:]
    #     # self.I_se=sol.y[5,:]
    #     # self.I_cr=sol.y[6,:]
    #     # self.H_in=sol.y[7,:]
    #     # self.H_cr=sol.y[8,:]
    #     # self.H_out=sol.y[9,:]
    #     # self.V=sol.y[10,:]
    #     # self.D=sol.y[11,:]
    #     # self.B=sol.y[12,:]
    #     # self.R=sol.y[13,:]

    #     return(sol)