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


class SEIRHVD:  
    """
        SEIRHVD Object:
        Construction:
            SEIRHVD(tsim,beta,mu,alpha,k=0,Htot=30,Vtot=20)
        with:
            tsim: Simulation time
            beta: Transmition rate
            mu: E/I initial rate
            alpha: Movility function (quarantine object function)
            k: Saturation Kynetics Dynamics factor
            Htot: Hospital capacity, either an int or a function(t)
            Vtot: VMI capacity, either an int or a function(t)
    """
    def __init__(self,tsim,beta,mu,alpha,k=0,Htot=30,Vtot=20,H0=0,V0=0,B0=0,D0=0,R0=0,I0=100,I_d0=10,I_ac0=100,SeroPrevFactor=1,expinfection=0,population=1000000,RealIC=None, initdate = None,Imi_det = 1,Ias_det = 1,SimIC=None):

        self.tsim = tsim 
        self.beta = beta 
        self.mu = mu 
        self.k = k
        self.alpha = alpha        

        self.SeroPrevFactor = SeroPrevFactor
        self.expinfection = expinfection

        self.Imi_det = Imi_det # Fraction of mild infected detected
        self.Ias_det = Ias_det # Fraction of asymptomatic infected detected

        self.SimIC = SimIC
        """
        I0 = Imi_det*Imi + Ias_det*Ias + Icr + Ise 
        Itot = Imi +  Icr + Ise  + Iac 
        Ias = Itot * pas
        Imi = Itot * pmi
        Ise = Itot * pse
        Icr = Itot * pcr

        I0 = Imi_det*Itot*pmi + Itot*(pcr + pse) + Ias_det*Itot*pas
        I0 = Itot(Imi_det*pmi + Ias_det*pas + pcr + pse)
        
        Itot = I0/(Imi_det*pmi + Ias_det*pas + pcr + pse)
        """
        self.gw=20

        self.initdate = initdate 
        # Initial Conditions:        
        #self.setinitvalues()  
        if RealIC:
            IC = RealIC
            # Use initial conditions object/dictionary with realworld imported initial conditions
            
            # Aproximate Hospitals capacity:
            Hcmodel = np.poly1d(np.polyfit(IC.sochimi_tr, IC.Hr_tot, 4))
            tsat = IC.sochimi_tr[-1]
            Hmax = np.mean(IC.Hr_tot[-10:])
            self.Htot=lambda t: Hcmodel(t)*(1-expit(t-tsat)) + expit(t-tsat)*Hmax  

            Vcmodel = np.poly1d(np.polyfit(IC.sochimi_tr, IC.Vr_tot, 4))
            tsat = IC.sochimi_tr[-1]
            Vmax = np.mean(IC.Vr_tot[-10:])
            self.Vtot=lambda t: Vcmodel(t)*(1-expit(t-tsat)) + expit(t-tsat)*Vmax

            # Set Initial values
            self.H0 = IC.Hr[0]
            self.V = IC.Vr[0]
            self.B = IC.Br[0]
            self.D = IC.Br[1]-IC.Br[0]
            self.R = 0
            self.I0 = IC.Ir[0]
            self.I_d0 = IC.I_d_r[0]
            self.I_ac0 = IC.I_ac_r[0]
            self.population = IC.population
            self.initdate = IC.initdate

            self.Ise_D_d = 0
            self.Icr_D_d = 0
            self.Hse_D_d = 0
            self.V_D_d = 0

            self.Ise_D_ac = 0
            self.Icr_D_ac = 0
            self.Hse_D_ac = 0
            self.V_D_ac = 0

            self.R_d = 0

            self.Hse_d = 0
            self.Hout_d = 0
            self.V_d = 0

            self.Hse_ac = 0
            self.Hout_ac = 0
            self.V_ac = 0
            
            self.setparams()
            self.Einit = False
            self.setrelationalvalues()
            self.setequations()
            print('InitialCondition Object Data')

        elif SimIC:
            self.SimICinitdate = SimIC.initdate
            self.population = SimIC.population
            
            # New Susceptible:
            self.S = SimIC.S[-1] + SimIC.population*(1-SimIC.SeroPrevFactor)*self.SeroPrevFactor
            self.N = SimIC.SeroPrevFactor*self.population + self.population*(1-SimIC.SeroPrevFactor)*self.SeroPrevFactor #Past simulation + added now

            # Exposed: 
            self.E = SimIC.E[-1]
            self.E_d = SimIC.E_d[-1]
            self.E_ac = SimIC.E_ac[-1]

            self.I = SimIC.I[-1]
            self.I_d = SimIC.I_d[-1]
            self.I_ac = SimIC.I_ac[-1]

            self.Ias = SimIC.Ias[-1]
            self.Imi = SimIC.Imi[-1]
            self.Ise = SimIC.Ise[-1]
            self.Icr = SimIC.Icr[-1]

            self.Ias_d = SimIC.Ias_d[-1]
            self.Imi_d = SimIC.Imi_d[-1]
            self.Ise_d = SimIC.Ise_d[-1]
            self.Icr_d = SimIC.Icr_d[-1]

            self.Ias_ac = SimIC.Ias_ac[-1]
            self.Imi_ac = SimIC.Imi_ac[-1]
            self.Ise_ac = SimIC.Ise_ac[-1]
            self.Icr_ac = SimIC.Icr_ac[-1]

            self.R = SimIC.R[-1]
            self.R_d = SimIC.R_d[-1]

            self.Hse = SimIC.Hse[-1]
            self.Hout = SimIC.Hout[-1]
            self.V = SimIC.V[-1]
            
            self.Hse_d = SimIC.Hse_d[-1]
            self.Hout_d = SimIC.Hout_d[-1]
            self.V_d = SimIC.V_d[-1]

            self.Hse_ac = SimIC.Hse_ac[-1]            
            self.Hout_ac = SimIC.Hout_ac[-1]            
            self.V_ac = SimIC.V_ac[-1]              

            self.R = SimIC.R[-1]
            self.R_d = SimIC.R_d[-1] 

            self.D = SimIC.D[-1]
            self.B = SimIC.B[-1]

            self.Ise_D_d = SimIC.Ise_D_d[-1]
            self.Icr_D_d = SimIC.Icr_D_d[-1]
            self.Hse_D_d = SimIC.Hse_D_d[-1]
            self.V_D_d = SimIC.V_D_d[-1]

            self.Ise_D_ac = SimIC.Ise_D_ac[-1]
            self.Icr_D_ac = SimIC.Icr_D_ac[-1]
            self.Hse_D_ac = SimIC.Hse_D_ac[-1]
            self.V_D_ac = SimIC.V_D_ac[-1]

            self.Einit = True

            # Falta trasladarlo en el tiempo
            self.T_delta = (self.initdate - SimIC.initdate).days

            self.Htot = self.Htot_SimIC(SimIC)
            self.Vtot = self.Vtot_SimIC(SimIC)

            self.alpha = self.alpha_SimIC(SimIC)

            self.setparams()
            self.setequations()              




        else:
            self.H0 = H0
            self.V = V0
            self.B = B0
            self.D = D0
            self.R = R0
            self.I0 = I0
            self.I_d0 = I_d0
            self.I_ac0 = I_ac0
            self.population = population

            # Build Hospital Capacities
            if type(Htot)==int:
                self.Htot = np.poly1d(Htot) 
            else:
                self.Htot = Htot

            if type(Vtot)==int:
                self.Vtot = np.poly1d(Vtot) 
            else:
                self.Vtot = Vtot            

            self.Ise_D_d = 0
            self.Icr_D_d = 0
            self.Hse_D_d = 0
            self.V_D_d = 0

            self.Ise_D_ac = 0
            self.Icr_D_ac = 0
            self.Hse_D_ac = 0
            self.V_D_ac = 0

            self.Hse_d = 0
            self.Hout_d = 0
            self.V_d = 0

            self.Hse_ac = 0
            self.Hout_ac = 0
            self.V_ac = 0            

            self.R_d = 0

            self.setparams()
            self.Einit = False
            self.setrelationalvalues()
            self.setequations() 


    def alpha_SimIC(self,SimIC):
        def funct(t):
            return SimIC.alpha(t+self.T_delta)
        return funct    


    def Htot_SimIC(self,SimIC):
        def funct(t):
            return SimIC.Htot(t+self.T_delta)
        return funct
        

    def Vtot_SimIC(self,SimIC):
        def funct(t):
            return SimIC.Vtot(t+self.T_delta)
        return funct
        #return(SimIC.Vtot(t-self.T_delta))               

    def setnewparams(self):
        self.setequations()
        if not self.SimIC:
            self.setrelationalvalues()
        print('Compartimental model State parameters changed')


    def setequations(self):
        """
        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #
        """
        # dVariable/dt = sum(prob_i/in_time_i*in_State_i,i in in states) - sum(prob_i/out_time_i*out_State_i,i in in states) 
        
        # --------------------------- #
        #        Susceptibles         #
        # --------------------------- #        
        
        # 0) dS/dt:
        self.dS=lambda t,S,E,Ias,Imi,Ise,Icr,R,D: -self.alpha(t)*self.beta*S*(self.expinfection*E+Ias+Imi+Ise+Icr)/(self.N+self.k*(Ias+Imi+Ise+Icr))-self.betaD*D+self.eta*R
        
        # --------------------------- #
        #           Exposed           #
        # --------------------------- #        
        
        # 1) dE_as/dt
        self.dE = lambda t,S,E,Ias,Imi,Ise,Icr: self.alpha(t)*self.beta*S*(self.expinfection*E+Ias+Imi+Ise+Icr)/(self.N+self.k*(Ias+Imi+Ise+Icr)) \
            -self.pE_Ias/self.tE_Ias*E -self.pE_Imi/self.tE_Imi*E-self.pE_Ise/self.tE_Ise*E-self.pE_Icr/self.tE_Icr*E
 
        # 2) Daily dE_as/dt
        self.dE_d = lambda t,S,E,E_d,Ias,Imi,Ise,Icr: self.alpha(t)*self.beta*S*(self.expinfection*E+Ias+Imi+Ise+Icr)/(self.N+self.k*(Ias+Imi+Ise+Icr)) - E_d
       
        # 3) Accumulated dE_as/dt
        self.dE_ac = lambda t,S,E,Ias,Imi,Ise,Icr: self.alpha(t)*self.beta*S*(self.expinfection*E+Ias+Imi+Ise+Icr)/(self.N+self.k*(Ias+Imi+Ise+Icr))


        # --------------------------- #
        #           Infected          #
        # --------------------------- #                
        
        #  --- Active --- #

        # 4) Asymptomatic dIas/dt
        self.dIas=lambda t,E,Ias: self.pE_Ias/self.tE_Ias*E-self.pIas_R/self.tIas_R*Ias
        # 5) Mild  dImi/dt
        self.dImi=lambda t,E,Imi: self.pE_Imi/self.tE_Imi*E-self.pImi_R/self.tImi_R*Imi
        # 6) Serious dIse/dt: Esy -  
        self.dIse=lambda t,E,Ise,Hse,Hout: self.pE_Ise/self.tE_Ise*E-self.pIse_Hse/self.tIse_Hse*Ise*(self.h_sat(Hse,Hout,t))\
            -self.pIse_D/self.tIse_D*Ise*(1-self.h_sat(Hse,Hout,t))            
        # 7) Critical  dIcr/dt
        self.dIcr=lambda t,E,Icr,V: self.pE_Icr/self.tE_Icr*E\
            -self.pIcr_V/self.tIcr_V*Icr*(self.v_sat(V,t)) - self.pIcr_D/self.tIcr_D*Icr*(1-self.v_sat(V,t))

        
        # 8-11) Daily Infected
        self.dIas_d = lambda t,E,Ias_d: self.pE_Ias/self.tE_Ias*E - Ias_d 
        self.dImi_d = lambda t,E,Imi_d: self.pE_Imi/self.tE_Imi*E - Imi_d 
        self.dIse_d = lambda t,E,Ise_d: self.pE_Ise/self.tE_Ise*E - Ise_d
        self.dIcr_d = lambda t,E,Icr_d: self.pE_Icr/self.tE_Icr*E - Icr_d
        #self.dI_d = lambda t,E,I_d: self.pE_Ias/self.tE_Ias*E + self.pE_Imi/self.tE_Imi*E + self.pE_Ise/self.tE_Ise*E + self.pE_Icr/self.tE_Icr*E - I_d

        # 12-15) Accummulated Infected
        self.dIas_ac = lambda t,E: self.pE_Ias/self.tE_Ias*E
        self.dImi_ac = lambda t,E: self.pE_Imi/self.tE_Imi*E
        self.dIse_ac = lambda t,E: self.pE_Ise/self.tE_Ise*E
        self.dIcr_ac = lambda t,E: self.pE_Icr/self.tE_Icr*E 
        #self.dI_ac = lambda t,E: self.pE_Ias/self.tE_Ias*E + self.pE_Imi/self.tE_Imi*E + self.pE_Ise/self.tE_Ise*E + self.pE_Icr/self.tE_Icr*E

        # --------------------------- #
        #        Hospitalized         #
        # --------------------------- #  
        
        # 16) dHse/dt: Serious Infected Hospitalized
        self.dHse=lambda t,Ise,Hse,Hout,V: self.pIse_Hse/self.tIse_Hse*Ise*(self.h_sat(Hse,Hout,t)) - self.pHse_V/self.tHse_V*Hse*(self.v_sat(V,t)) \
             - self.pHse_D/self.tHse_D*Hse*(1- self.v_sat(V,t)) - self.pHse_R/self.tHse_R*Hse
            

        # 17) dHout/dt: Hospitalized Recovering after VMI
        self.dHout=lambda t,Hout,V: self.pV_Hout/self.tV_Hout*V-self.pHout_R/self.tHout_R*Hout

        # 18) dV/dt: Ventilator use
        self.dV=lambda t,Icr,Hse,V: self.pIcr_V/self.tIcr_V*Icr*(self.v_sat(V,t))+ self.pHse_V/self.tHse_V*Hse*(self.v_sat(V,t)) \
            -  self.pV_Hout/self.tV_Hout*V - self.pV_D/self.tV_D*V

        # Daily 

        # 19) dHse/dt: Serious Infected Hospitalized
        self.dHse_d=lambda t,Ise,Hse,Hout,Hse_d: self.pIse_Hse/self.tIse_Hse*Ise*(self.h_sat(Hse,Hout,t)) - Hse_d
            

        # 20) dHout/dt: Hospitalized Recovering after VMI
        self.dHout_d=lambda t,V,Hout_d: self.pV_Hout/self.tV_Hout*V - Hout_d

        # 21) dV/dt: Ventilator use
        self.dV_d=lambda t,Icr,Hse,V,V_d: self.pIcr_V/self.tIcr_V*Icr*(self.v_sat(V,t))+ self.pHse_V/self.tHse_V*Hse*(self.v_sat(V,t)) - V_d
            

        # Accumulated:
        # 22) dHse/dt: Serious Infected Hospitalized
        self.dHse_ac=lambda t,Ise,Hse,Hout: self.pIse_Hse/self.tIse_Hse*Ise*(self.h_sat(Hse,Hout,t))             

        # 23) dHout/dt: Hospitalized Recovering after VMI
        self.dHout_ac=lambda t,V: self.pV_Hout/self.tV_Hout*V 

        # 24) dV/dt: Ventilator use
        self.dV_ac=lambda t,Icr,Hse,V: self.pIcr_V/self.tIcr_V*Icr*(self.v_sat(V,t))+ self.pHse_V/self.tHse_V*Hse*(self.v_sat(V,t))


        # --------------------------- #
        #         Recovered           #
        # --------------------------- #  
        
        # 25) dR/dt
        self.dR=lambda t,Ias,Imi,Hse,Hout,R: self.pIas_R/self.tIas_R*Ias + self.pImi_R/self.tImi_R*Imi + \
            self.pHout_R/self.tHout_R*Hout + self.pHse_R/self.tHse_R*Hse - self.eta*R

        # 26) Daily recovered
        self.dR_d=lambda t,Ias,Imi,Hse,Hout,R,R_d: self.pIas_R/self.tIas_R*Ias + self.pImi_R/self.tImi_R*Imi + \
            self.pHout_R/self.tHout_R*Hout + self.pHse_R/self.tHse_R*Hse - R_d


        # --------------------------- #
        #           Deaths            #
        # --------------------------- #         
        
        # 27) dD/dt: Death Rate
        self.dD=lambda t,Ise,Icr,Hse,Hout,V,D: self.pIse_D/self.tIse_D*Ise*(1-self.h_sat(Hse,Hout,t)) + self.pIcr_D/self.tIcr_D*Icr*(1-self.v_sat(V,t)) + \
            self.pHse_D/self.tHse_D*Hse*(1- self.v_sat(V,t)) + self.pV_D/self.tV_D*V - self.pD_B/self.tD_B*D

        # 28) dB/dt: Bury rate
        self.dB=lambda t,D: self.pD_B/self.tD_B*D
        
        # 29-32) Daily Deads         
        self.dIse_D_d = lambda t,Ise,Hse,Hout,Ise_D_d: self.pIse_D/self.tIse_D*Ise*(1-self.h_sat(Hse,Hout,t)) - Ise_D_d
        self.dIcr_D_d = lambda t,Icr,Hse,Hout,V,Icr_D_d:  self.pIcr_D/self.tIcr_D*Icr*(1-self.v_sat(V,t)) - Icr_D_d
        self.dHse_D_d = lambda t,Hse,V,Hse_D_d: self.pHse_D/self.tHse_D*Hse*(1- self.v_sat(V,t)) - Hse_D_d
        self.dV_D_d   = lambda t,V,V_D_d: self.pV_D/self.tV_D*V - V_D_d
        
        # 33-36) Accumulated Deads ṕer cause
        self.dIse_D_ac = lambda t,Ise,Hse,Hout: self.pIse_D/self.tIse_D*Ise*(1-self.h_sat(Hse,Hout,t)) 
        self.dIcr_D_ac = lambda t,Icr,Hse,Hout,V:  self.pIcr_D/self.tIcr_D*Icr*(1-self.v_sat(V,t)) 
        self.dHse_D_ac = lambda t,Hse,V: self.pHse_D/self.tHse_D*Hse*(1- self.v_sat(V,t)) 
        self.dV_D_ac   = lambda t,V: self.pV_D/self.tV_D*V 


    # UCI and UTI beds saturation function
    def h_sat(self,Hse,Hout,t):
        return(expit(-self.gw*(Hse+Hout-self.Htot(t))))
    # Ventilators Saturation Function    
    def v_sat(self,V,t):
        return(expit(-self.gw*(V-self.Vtot(t))))

    def setparams(self):
        self.pE_Ias = 0.4  # Transition from exposed to Asymptomatic Infected
        self.tE_Ias = 5.0

        self.pE_Imi = 0.55 # Transition from exposed to  Mild Infected
        self.tE_Imi = 5.0

        self.pE_Icr = 0.01666 # Transition from exposed to  Critical Infected
        self.tE_Icr = 3.0

        self.pE_Ise = 0.03334 ## Transition from exposed to  Serious Infected
        self.tE_Ise = 3.0

        self.pIas_R = 1.0   # Transition from Asymptomatic Infected to Recovered
        self.tIas_R = 10.0 

        self.pImi_R = 1.0  # Transition from Mild Infected to Recovered
        self.tImi_R = 15.0

        self.pIse_Hse = 1.0 # Transition from Serious Infected to Serious Hospitalized (When Hospital capacity is not saturated)
        self.tIse_Hse = 3.0 

        self.pIse_D = 1.0  # Transition from Serious Infected to Death (When Hospital capacity is saturated)
        self.tIse_D = 3.0         

        self.pIcr_V = 1.0  # Transition from Critical Infected to Ventilator (When Ventilators capacity is not saturated)
        self.tIcr_V = 3.0 

        self.pIcr_D = 1.0  # Transition from Serious Infected to Death (When Ventilators capacity is saturated)
        self.tIcr_D = 3.0         

        self.pHse_R = 0.97 # Transition from Serious Hospitalized to Recovered
        self.tHse_R = 11.0

        self.pHse_V = 0.03 # Transition from Serious Hospitalized to Ventilators (When Ventilators capacity is not saturated)
        self.tHse_V = 3.0

        self.pHse_D = 0.03 # Transition from Serious Hospitalized to Death (When Ventilators capacity is saturated)
        self.tHse_D = 3.0        

        self.pV_Hout = 0.5  # Transition from Ventilators to Hospital Recovery (Hout) 
        self.tV_Hout = 15.0

        self.pV_D = 0.5 # Transition from Ventilators to Death
        self.tV_D = 15.0

        self.pHout_R = 1.0 # Transition from Hospital Recovery (Hout) to Recovered
        self.tHout_R = 4.0

        self.pD_B = 1.0 # Transition from Dead to buried
        self.tD_B = 1.0 

        self.betaD = 0 # Contagion by deads rate
        self.eta = 0.0 # Immunity loss rate


        # ------------------- #
        #  Valores Iniciales  #
        # ------------------- #
            
    def setinitvalues(self):
        # 15 de Mayo
        self.I_act0 = 12642

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
        self.Ias_ac = 0
        self.Imi_ac = 0
        self.Ise_ac = 0
        self.Icr_ac = 0

        # Deaths
        self.H_crD = 0
        self.VD = 0
        self.Ise_D = 0
        self.IcrD = 0

        # Daily Infected
        self.Ias_d = 0
        self.Imi_d = 0
        self.Ise_d = 0
        self.Icr_d = 0

        # Daily Deaths
        self.H_crD_d = 0
        self.VD_d = 0
        self.Ise_D_d = 0
        self.IcrD_d = 0

        # Initial Infected proportion
        self.Ias_prop = 0.35
        self.Imi_prop = 0.63
        self.Icr_prop = 0.007
        self.Ise_prop = 0.013

        self.setrelationalvalues()

    def setrelationalvalues(self):
        self.I = self.I0/(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr )
        # Active infected
        self.Ias= self.pE_Ias*self.I
        self.Imi= self.pE_Imi*self.I
        self.Icr= self.pE_Icr*self.I
        self.Ise = self.pE_Ise*self.I

        self.I_d = self.I_d0/(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr )
        self.Ias_d = self.pE_Ias*self.I_d
        self.Imi_d = self.pE_Imi*self.I_d
        self.Icr_d = self.pE_Icr*self.I_d
        self.Ise_d = self.pE_Ise*self.I_d

        self.I_ac = self.I_ac0/(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr)
        self.Ias_ac = self.pE_Ias*self.I_ac
        self.Imi_ac = self.pE_Imi*self.I_ac
        self.Icr_ac = self.pE_Icr*self.I_ac
        self.Ise_ac = self.pE_Ise*self.I_ac
        
        # Exposed
        if not self.Einit:
            self.E = self.mu*self.I
            self.E_d=self.mu*(self.I_d)                
            self.E_ac=self.mu*(self.I_ac)
        # Hospitalizados        
        self.Hse = self.H0*self.pE_Ise/(self.pE_Ise+self.pE_Icr)
        self.Hout = self.H0*self.pE_Icr/(self.pE_Ise+self.pE_Icr)
       
        # Valores globales
        self.N = self.SeroPrevFactor*self.population
        self.S = self.N-self.H0-self.V-self.D-self.E-(self.Ias+self.Icr+self.Ise+self.Imi)        
        #self.I = self.I

        #constructor of SEIR class elements, it's initialized when a parameter
        #miminization is performed to adjust the best setting of the actual infected

    """
    def calculateindicators(self):
        self.R_ef
        self.SHFR
        # Peak
        self.peak
        self.peakindex
        self.peak_t
        self.peak_date
        # Saturation dates
        self.Hsat_t
        self.Hsat_date
        self.VMIsat_t
        self.VMIsat_date

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
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.
        print('Import odeint')
        from scikits.odes.odeint import odeint


        if(not isinstance(self.S, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                E0=self.mu*(self.I)
                E_d0=self.mu*(self.I_d)
                E_ac0=self.mu*(self.I_ac)
            else:
                E0 = self.E
                E_d0 = self.mu*(self.I_d)
                E_ac0 = self.mu*(self.I_ac)                

            S0=self.S

            E0=self.E
            E_d0=self.E_d
            E_ac0=self.E_ac

            Ias0=self.Ias
            Imi0=self.Imi
            Ise0=self.Ise
            Icr0=self.Icr

            Ias_d0=self.Ias_d
            Imi_d0=self.Imi_d
            Ise_d0=self.Ise_d
            Icr_d0=self.Icr_d   

            Ias_ac0=self.Ias_ac
            Imi_ac0=self.Imi_ac
            Ise_ac0=self.Ise_ac
            Icr_ac0=self.Icr_ac
        
            Hse0=self.Hse
            Hout0=self.Hout
            V0=self.V

            Hse_d0= self.Hse_d
            Hout_d0= self.Hout_d
            V_d0= self.V_d

            Hse_ac0= self.Hse_ac
            Hout_ac0= self.Hout_ac
            V_ac0= self.V_ac

            R0=self.R
            R_d0=self.R_d

            D0=self.D
            B0=self.B

            Ise_D_d0 = self.Ise_D_d
            Icr_D_d0 = self.Icr_D_d
            Hse_D_d0 = self.Hse_D_d
            V_D_d0 = self.V_D_d

            Ise_D_ac0 = self.Ise_D_ac
            Icr_D_ac0 = self.Icr_D_ac
            Hse_D_ac0 = self.Hse_D_ac
            V_D_ac0 = self.V_D_ac



            self.t=np.arange(t0,T+h,h)
            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            E0 = self.E
            E_d0 = self.E_d
            E_ac0 ==self.E_ac

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


            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)

        else:
            return()

        
        def model_SEIR_graph(t,y,ydot):
            
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3],y[4],y[5],y[19],y[20])

            ydot[1]=self.dE(t,y[0],y[1],y[2],y[3],y[4],y[5])
            ydot[2]=self.dE_d(t,y[0],y[1],y[2],y[4],y[5],y[6],y[7])
            ydot[3]=self.dE_ac(t,y[0],y[1],y[4],y[5],y[6],y[7])
                                
            ydot[4]=self.dIas(t,y[1],y[4])
            ydot[5]=self.dImi(t,y[1],y[5])
            ydot[6]=self.dIse(t,y[1],y[6],y[16],y[17])
            ydot[7]=self.dIcr(t,y[1],y[7],y[18])

            ydot[8]=self.dIas_d(t,y[1],y[8])
            ydot[9]=self.dImi_d(t,y[1],y[9]) 
            ydot[10]=self.dIse_d(t,y[1],y[10]) 
            ydot[11]=self.dIcr_d(t,y[1],y[11])            

            ydot[12]=self.dIas_ac(t,y[1])
            ydot[13]=self.dImi_ac(t,y[1]) 
            ydot[14]=self.dIse_ac(t,y[1]) 
            ydot[15]=self.dIcr_ac(t,y[1])
            
            ydot[16]=self.dHse(t,y[6],y[16],y[17],y[18])            
            ydot[17]=self.dHout(t,y[17],y[18])
            ydot[18]=self.dV(t,y[7],y[16],y[18])

            ydot[19]=self.dHse_d(t,y[6],y[16],y[17],y[19])            
            ydot[20]=self.dHout_d(t,y[18],y[20])
            ydot[21]=self.dV_d(t,y[7],y[16],y[18],y[21])

            ydot[22]=self.dHse_ac(t,y[6],y[16],y[17])
            ydot[23]=self.dHout_ac(t,y[18])
            ydot[24]=self.dV_ac(t,y[7],y[16],y[18])

            ydot[25]=self.dR(t,y[4],y[5],y[16],y[17],y[25])
            ydot[26]=self.dR_d(t,y[4],y[5],y[16],y[17],y[25],y[26])

            ydot[27]=self.dD(t,y[6],y[7],y[16],y[17],y[18],y[27])
            ydot[28]=self.dB(t,y[27])

            ydot[29]=self.dIse_D_d(t,y[6],y[16],y[17],y[29])
            ydot[30]=self.dIcr_D_d(t,y[7],y[16],y[17],y[18],y[30])
            ydot[31]=self.dHse_D_d(t,y[16],y[18],y[31])
            ydot[32]=self.dV_D_d(t,y[18],y[32])

            ydot[33]=self.dIse_D_ac(t,y[6],y[16],y[17])
            ydot[34]=self.dIcr_D_ac(t,y[7],y[16],y[17],y[18])
            ydot[35]=self.dHse_D_ac(t,y[16],y[18])
            ydot[36]=self.dV_D_ac(t,y[18])



        initcond = np.array([S0,E0,E_d0,E_ac0,Ias0,Imi0,Ise0,Icr0,Ias_d0,Imi_d0,Ise_d0,Icr_d0,Ias_ac0,Imi_ac0,Ise_ac0,Icr_ac0,Hse0,Hout0,V0,Hse_d0,Hout_d0,V_d0,Hse_ac0,
            Hout_ac0,V_ac0,R0,R_d0,D0,B0,Ise_D_d0,Icr_D_d0,Hse_D_d0,V_D_d0,Ise_D_ac0,Icr_D_ac0,Hse_D_ac0,V_D_ac0])
            
        print('Solving ODE')
        sol = odeint(model_SEIR_graph, self.t, initcond,method='admo')
        
        self.t=sol.values.t 
        
        self.S=sol.values.y[:,0]

        self.E=sol.values.y[:,1]
        self.E_d=sol.values.y[:,2]
        self.E_ac=sol.values.y[:,3]

        self.Ias=sol.values.y[:,4]
        self.Imi=sol.values.y[:,5]
        self.Ise=sol.values.y[:,6]
        self.Icr=sol.values.y[:,7]

        self.Ias_d=sol.values.y[:,8]
        self.Imi_d=sol.values.y[:,9]
        self.Ise_d=sol.values.y[:,10]
        self.Icr_d=sol.values.y[:,11]

        self.Ias_ac=sol.values.y[:,12]
        self.Imi_ac=sol.values.y[:,13]
        self.Ise_ac=sol.values.y[:,14]
        self.Icr_ac=sol.values.y[:,15]

        self.Hse=sol.values.y[:,16]
        self.Hout=sol.values.y[:,17]
        self.V=sol.values.y[:,18]

        self.Hse_d=sol.values.y[:,19]
        self.Hout_d=sol.values.y[:,20]
        self.V_d=sol.values.y[:,21]

        self.Hse_ac=sol.values.y[:,22]
        self.Hout_ac=sol.values.y[:,23]
        self.V_ac=sol.values.y[:,24]                

        self.R=sol.values.y[:,25]
        self.R_d=sol.values.y[:,26]

        self.D=sol.values.y[:,27]
        self.B=sol.values.y[:,28]
        
        self.Ise_D_d=sol.values.y[:,29]
        self.Icr_D_d=sol.values.y[:,30]
        self.Hse_D_d=sol.values.y[:,31]
        self.V_D_d=sol.values.y[:,32]

        self.Ise_D_ac=sol.values.y[:,33]
        self.Icr_D_ac=sol.values.y[:,34]
        self.Hse_D_ac=sol.values.y[:,35]
        self.V_D_ac=sol.values.y[:,36]



        self.I = self.Ias + self.Imi + self.Ise + self.Icr
        self.I_d = self.Ias_d + self.Imi_d + self.Ise_d + self.Icr_d
        self.I_ac = self.Ias_ac + self.Imi_ac + self.Ise_ac + self.Icr_ac

        self.H = self.Hse + self.Hout
        self.H_d = self.Hse_d + self.Hout_d
        self.H_ac = self.Hse_ac + self.Hout_ac

        self.H_sat = [self.h_sat(self.Hse[i],self.Hout[i],self.t[i]) for i in range(len(self.t))]
        self.V_sat = [self.v_sat(self.V[i],self.t[i]) for i in range(len(self.t))]

        self.V_cap = [self.Vtot(i) for i in self.t]
        self.H_cap = [self.Htot(i) for i in self.t]

        #Cálculo de la fecha del Peak  
        self.peakindex = np.where(self.I==max(self.I))[0][0]
        self.peak = max(self.I)
        self.peak_t = self.t[self.peakindex]
        if self.initdate:
            self.peak_date = self.initdate+timedelta(days=round(self.peak_t))

        # Detected Cases
        self.I_det = self.I*(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr )
        self.I_d_det = self.I_d*(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr )
        self.I_ac_det = self.I_ac*(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr )

        return(sol)

    def integr_sci(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.

        if(not isinstance(self.S, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                E0=self.mu*(self.I)
                E_d0=self.mu*(self.I_d)
                E_ac0=self.mu*(self.I_ac)
            else:
                E0 = self.E
                E_d0 = self.mu*(self.I_d)
                E_ac0 = self.mu*(self.I_ac)

            S0=self.S

            E0=self.E
            E_d0=self.E_d
            E_ac0=self.E_ac

            Ias0=self.Ias
            Imi0=self.Imi
            Ise0=self.Ise
            Icr0=self.Icr

            Ias_d0=self.Ias_d
            Imi_d0=self.Imi_d
            Ise_d0=self.Ise_d
            Icr_d0=self.Icr_d   

            Ias_ac0=self.Ias_ac
            Imi_ac0=self.Imi_ac
            Ise_ac0=self.Ise_ac
            Icr_ac0=self.Icr_ac
        
            Hse0=self.Hse
            Hout0=self.Hout
            V0=self.V

            Hse_d0= self.Hse_d
            Hout_d0= self.Hout_d
            V_d0= self.V_d

            Hse_ac0= self.Hse_ac
            Hout_ac0= self.Hout_ac
            V_ac0= self.V_ac

            R0=self.R
            R_d0=self.R_d

            D0=self.D
            B0=self.B

            Ise_D_d0 = self.Ise_D_d
            Icr_D_d0 = self.Icr_D_d
            Hse_D_d0 = self.Hse_D_d
            V_D_d0 = self.V_D_d

            Ise_D_ac0 = self.Ise_D_ac
            Icr_D_ac0 = self.Icr_D_ac
            Hse_D_ac0 = self.Hse_D_ac
            V_D_ac0 = self.V_D_ac

            self.t=np.arange(t0,T+h,h)
            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            E0 = self.E
            E_d0 = self.E_d
            E_ac0 ==self.E_ac

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

            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)


        else:
            return()
            

        
        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dS(t,y[0],y[1],y[2],y[3],y[4],y[5],y[19],y[20])

            ydot[1]=self.dE(t,y[0],y[1],y[2],y[3],y[4],y[5])
            ydot[2]=self.dE_d(t,y[0],y[1],y[2],y[4],y[5],y[6],y[7])
            ydot[3]=self.dE_ac(t,y[0],y[1],y[4],y[5],y[6],y[7])
                                
            ydot[4]=self.dIas(t,y[1],y[4])
            ydot[5]=self.dImi(t,y[1],y[5])
            ydot[6]=self.dIse(t,y[1],y[6],y[16],y[17])
            ydot[7]=self.dIcr(t,y[1],y[7],y[18])

            ydot[8]=self.dIas_d(t,y[1],y[8])
            ydot[9]=self.dImi_d(t,y[1],y[9]) 
            ydot[10]=self.dIse_d(t,y[1],y[10]) 
            ydot[11]=self.dIcr_d(t,y[1],y[11])            

            ydot[12]=self.dIas_ac(t,y[1])
            ydot[13]=self.dImi_ac(t,y[1]) 
            ydot[14]=self.dIse_ac(t,y[1]) 
            ydot[15]=self.dIcr_ac(t,y[1])
            
            ydot[16]=self.dHse(t,y[6],y[16],y[17],y[18])            
            ydot[17]=self.dHout(t,y[17],y[18])
            ydot[18]=self.dV(t,y[7],y[16],y[18])

            ydot[19]=self.dHse_d(t,y[6],y[16],y[17],y[19])            
            ydot[20]=self.dHout_d(t,y[18],y[20])
            ydot[21]=self.dV_d(t,y[7],y[16],y[18],y[21])

            ydot[22]=self.dHse_ac(t,y[6],y[16],y[17])
            ydot[23]=self.dHout_ac(t,y[18])
            ydot[24]=self.dV_ac(t,y[7],y[16],y[18])

            ydot[25]=self.dR(t,y[4],y[5],y[16],y[17],y[25])
            ydot[26]=self.dR_d(t,y[4],y[5],y[16],y[17],y[25],y[26])

            ydot[27]=self.dD(t,y[6],y[7],y[16],y[17],y[18],y[27])
            ydot[28]=self.dB(t,y[27])

            ydot[29]=self.dIse_D_d(t,y[6],y[16],y[17],y[29])
            ydot[30]=self.dIcr_D_d(t,y[7],y[16],y[17],y[18],y[30])
            ydot[31]=self.dHse_D_d(t,y[16],y[18],y[31])
            ydot[32]=self.dV_D_d(t,y[18],y[32])

            ydot[33]=self.dIse_D_ac(t,y[6],y[16],y[17])
            ydot[34]=self.dIcr_D_ac(t,y[7],y[16],y[17],y[18])
            ydot[35]=self.dHse_D_ac(t,y[16],y[18])
            ydot[36]=self.dV_D_ac(t,y[18])
                           
                                          
            return(ydot)

        initcond = np.array([S0,E0,E_d0,E_ac0,Ias0,Imi0,Ise0,Icr0,Ias_d0,Imi_d0,Ise_d0,Icr_d0,Ias_ac0,Imi_ac0,Ise_ac0,Icr_ac0,Hse0,Hout0,V0,Hse_d0,Hout_d0,V_d0,Hse_ac0,
            Hout_ac0,V_ac0,R0,R_d0,D0,B0,Ise_D_d0,Icr_D_d0,Hse_D_d0,V_D_d0,Ise_D_ac0,Icr_D_ac0,Hse_D_ac0,V_D_ac0])

        
        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA')
        
        self.t=sol.t 
        
        self.S=sol.y[0,:]
        self.E=sol.y[1,:]
        self.E_d=sol.y[2,:]
        self.E_ac=sol.y[3,:]
        self.Ias=sol.y[4,:]
        self.Imi=sol.y[5,:]
        self.Ise=sol.y[6,:]
        self.Icr=sol.y[7,:]
        self.Ias_d=sol.y[8,:]
        self.Imi_d=sol.y[9,:]
        self.Ise_d=sol.y[10,:]
        self.Icr_d=sol.y[11,:]
        self.Ias_ac=sol.y[12,:]
        self.Imi_ac=sol.y[13,:]
        self.Ise_ac=sol.y[14,:]
        self.Icr_ac=sol.y[15,:]
        self.Hse=sol.y[16,:]
        self.Hout=sol.y[17,:]
        self.V=sol.y[18,:]
        self.Hse_d=sol.y[19,:]
        self.Hout_d=sol.y[20,:]
        self.V_d=sol.y[21,:]
        self.Hse_ac=sol.y[22,:]
        self.Hout_ac=sol.y[23,:]
        self.V_ac=sol.y[24,:]
        self.R=sol.y[25,:]
        self.R_d=sol.y[26,:]
        self.D=sol.y[27,:]
        self.B=sol.y[28,:]
        self.Ise_D_d=sol.y[29,:]
        self.Icr_D_d=sol.y[30,:]
        self.Hse_D_d=sol.y[31,:]
        self.V_D_d=sol.y[32,:]
        self.Ise_D_ac=sol.y[33,:]
        self.Icr_D_ac=sol.y[34,:]
        self.Hse_D_ac=sol.y[35,:]
        self.V_D_ac=sol.y[36,:]


        self.I = self.Ias + self.Imi + self.Ise + self.Icr
        self.I_d = self.Ias_d + self.Imi_d + self.Ise_d + self.Icr_d
        self.I_ac = self.Ias_ac + self.Imi_ac + self.Ise_ac + self.Icr_ac

        self.H = self.Hse + self.Hout
        self.H_d = self.Hse_d + self.Hout_d
        self.H_ac = self.Hse_ac + self.Hout_ac

        self.H_sat = [self.h_sat(self.Hse[i],self.Hout[i],self.t[i]) for i in range(len(self.t))]
        self.V_sat = [self.v_sat(self.V[i],self.t[i]) for i in range(len(self.t))]

        self.V_cap = [self.Vtot(i) for i in self.t]
        self.H_cap = [self.Htot(i) for i in self.t]

        #Cálculo de la fecha del Peak  
        self.peakindex = np.where(self.I==max(self.I))[0][0]
        self.peak = max(self.I)
        self.peak_t = self.t[self.peakindex]
        if self.initdate:
            self.peak_date = self.initdate+timedelta(days=round(self.peak_t))

        # Detected Cases
        self.I_det = self.I*(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr )
        self.I_d_det = self.I_d*(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr )
        self.I_ac_det = self.I_ac*(self.Ias_det*self.pE_Ias + self.Imi_det*self.pE_Imi + self.pE_Ise + self.pE_Icr )

        return(sol)






