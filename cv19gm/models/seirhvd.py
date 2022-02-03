#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIRHVD Model
"""

import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd
import toml
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



class SEIRHVD:  
    """
        SEIR model object:
        Construction:
            SEIR(self, config = None, inputdata=None)

    """
    def __init__(self, config, inputdata=None,verbose = False,**kwargs):        

        if not config:
            print('Missing configuration file ')
            return None

        # ------------------------------- #
        #        Parameters Load          #
        # ------------------------------- #
        self.config = config  
        if verbose:
            print('Loading configuration file')
        cv19files.loadconfig(self,config,inputdata,**kwargs) # Load configuration file
        
        # Hospital capacity:        
        if verbose:
            print('Building hospital capacity model')        
        
        # Use data
        if type(self.cfg['initialconditions']['H_cap'])==str:            
            self.H_cap = cv19functions.polyfit(self.data[self.H_cap],time=self.data[self.H_cap+'_t'])            

        self.H_sat = cv19functions.saturation(self.H_cap)
        
        if verbose:
            print('Initializing parameters and variables')
        self.set_relational_values()
        if verbose:
            print('Building equations')            
        self.set_equations()
        
        self.solved = False

        if verbose:
            print('SEIRHVD object created')

    # ------------------- #
    #  Valores Iniciales  #
    # ------------------- #
    # Init Ev
   
    def set_relational_values(self):

        #Initial Population
        self.N0 = self.popfraction*self.population
           
        # Vaccinated Infected
        vacprop = (1-self.vac_eff(0))*self.Sv/(self.N0 - self.Sv - self.E - self.I - self.H - self.D - self.R )

        if not self.Iv:
            self.Iv = self.I*vacprop
            self.Iv_d = self.I_d*vacprop
            self.Iv_ac = self.I_ac*vacprop

         # Infected
        self.Icr = (self.I-self.Iv)*self.pE_Icr(0)
        self.Im = (self.I-self.Iv)*self.pE_Im(0)
        self.Icr_d = (self.I_d-self.Iv_d)*self.pE_Icr(0)
        self.Im_d = (self.I_d-self.Iv_d)*self.pE_Im(0)
        self.Icr_ac = (self.I_ac-self.Iv_ac)*self.pE_Icr(0)
        self.Im_ac = (self.I_ac-self.Iv_ac)*self.pE_Im(0)

        
        # Exposed
        #if not self.Einit:
        self.E = self.mu*self.I
        self.E_d=self.mu*self.I_d                
        self.E_ac=self.mu*self.I_ac
       
        # Exposed vaccinated
        #if not self.Einit:
        self.Ev = self.mu*self.Iv
        self.Ev_d=self.mu*self.Iv_d                
        self.Ev_ac=self.mu*self.Iv_ac       
        
        # Initial susceptible population
        self.S = self.N0 - self.Sv - self.E - self.Ev - self.I - self.Iv - self.H - self.D - self.R

                    

    def set_equations(self):
        """
        # --------------------------- #
        #    Diferential Ecuations    #
        # --------------------------- #
        Variables: 
        S: Susceptibles
        Sv: Vaccinated Susceptibles
        E: Exposed
        I_m: Asymptomatic + mild + severe infected
        I_cr: Critical Infected
        I_v Vaccinated Infected
        Phi: Integrated external flux
        S_f: Susceptible external Flux
        """

        # --------------------------- #
        #        Susceptibles         #
        # --------------------------- # 
       
        # 0) dS/dt:
        self.dS=lambda t,S,Im,Icr,Iv,R,Phi: - self.alpha(t)*S*(self.beta(t)*(Im+Icr)+self.beta_v(t)*Iv)/(self.N0+Phi) + self.pR_S(t)/self.tR_S(t)*R - self.vac_d(t) + self.S_f(t)
        
        # 1) dS_v/dt:
        self.dSv=lambda t,Sv,Im,Icr,Iv,R,Phi: -(1-self.vac_eff(t))*self.alpha(t)*Sv*(self.beta(t)*(Im+Icr)+self.beta_v(t)*Iv)/(self.N0+Phi) + self.vac_d(t) + self.Sv_f(t)
                
        # --------------------------- #
        #           Exposed           #
        # --------------------------- #        
        
        # 2) dE/dt
        self.dE = lambda t,S,E,Im,Icr,Iv,Phi: self.alpha(t)*S*(self.beta(t)*(Im+Icr)+self.beta_v(t)*Iv)/(self.N0+Phi) - E*(self.pE_Im(t)/self.tE_Im(t)+self.pE_Icr(t)/self.tE_Icr(t)) + self.E_f(t)
        
        # 3) dE_d/dt*-
        self.dE_d = lambda t,S,E_d,Im,Icr,Iv,Phi: self.alpha(t)*S*(self.beta(t)*(Im+Icr)+self.beta_v(t)*Iv)/(self.N0+Phi) - E_d

        # 4) dEv/dt
        self.dEv = lambda t,Sv,Ev,Im,Icr,Iv,Phi: (1-self.vac_eff(t))*self.alpha(t)*Sv*(self.beta(t)*(Im+Icr)+self.beta_v(t)*Iv)/(self.N0+Phi) - Ev/self.tEv_Iv(t) + self.Ev_f(t) 

        # 5) dEv_d/dt
        self.dEv_d = lambda t,Sv,Ev_d,Im,Icr,Iv,Phi: (1-self.vac_eff(t))*self.alpha(t)*Sv*(self.beta(t)*(Im+Icr)+self.beta_v(t)*Iv)/(self.N0+Phi) - Ev_d
        
        # --------------------------- #
        #           Infected          #
        # --------------------------- #                
        
        # 6) Mild Infected: Actuve
        self.dIm=lambda t,E,Im: self.pE_Im(t)/self.tE_Im(t)*E - 1/self.tIm_R(t)*Im + self.Im_f(t)

        # 7) Mild Infected: Daily
        self.dIm_d=lambda t,E,Im_d: self.pE_Im(t)/self.tE_Im(t)*E - Im_d

        # 8) Critical Infected: Active
        self.dIcr=lambda t,E,Icr: self.pE_Icr(t)/self.tE_Icr(t)*E - 1/self.tIcr_H(t)*Icr + self.Icr_f(t)

        # 9) Critical Infected: Daily
        self.dIcr_d=lambda t,E,Icr_d: self.pE_Icr(t)/self.tE_Icr(t)*E - Icr_d

        # 10) Vaccinated Infected: Active
        self.dIv=lambda t,Ev,Iv: Ev/self.tEv_Iv(t) - self.pIv_R(t)/self.tIv_R(t)*Iv - self.pIv_H(t)/self.tIv_H(t)*Iv + self.Iv_f(t)
        
        # 11) Vaccinated Infected: Daily
        self.dIv_d = lambda t,Ev,Iv_d: Ev/self.tEv_Iv(t) - Iv_d 

        # ---------------------------- #
        #        Hospitalized          #
        # ---------------------------- #  
        # 12) Hospitalized
        self.dH=lambda t,Icr,Iv,H: (1-self.H_sat(t,H))*(1/self.tIcr_H(t)*Icr + self.pIv_H(t)/self.tIv_H(t)*Iv) - self.pH_R(t)/self.tH_R(t)*H - self.pH_D(t)/self.tH_D(t)*H  + self.H_f(t)

        # 13) Hospitalized: Daily 
        self.dH_d=lambda t,Icr,Iv,H,H_d: (1-self.H_sat(t,H))*(1/self.tIcr_H(t)*Icr + self.pIv_H(t)/self.tIv_H(t)*Iv) - H_d

        # ------------------------- #
        #           Deaths          #
        # ------------------------- #  
        # 14) Deaths: total
        self.dD=lambda t,Icr,Iv,H: self.H_sat(t,H)*(1/self.tIcr_H(t)*Icr + self.pIv_H(t)/self.tIv_H(t)*Iv) + self.pH_D(t)/self.tH_D(t)*H  + self.D_f(t)

        # 15) Deaths: Daily
        self.dD_d=lambda t,Icr,Iv,H,D_d: self.H_sat(t,H)*(1/self.tIcr_H(t)*Icr + self.pIv_H(t)/self.tIv_H(t)*Iv) + self.pH_D(t)/self.tH_D(t)*H  - D_d

        # --------------------------- #
        #         Recovered           #
        # --------------------------- #  
        
        # 16) Total recovered
        self.dR=lambda t,Im,Iv,H,R: self.pIv_R(t)/self.tIv_R(t)*Iv + 1/self.tIm_R(t)*Im + self.pH_R(t)/self.tH_R(t)*H - self.pR_S(t)/self.tR_S(t)*R + self.R_f(t)

        # 17) Total recovered
        self.dR_d=lambda t,Im,Iv,H,R_d: self.pIv_R(t)/self.tIv_R(t)*Iv + 1/self.tIm_R(t)*Im + self.pH_R(t)/self.tH_R(t)*H - R_d

        # 18) External Flux:
        self.dPhi = lambda t: self.S_f(t) + self.Sv_f(t) + self.E_f(t) + self.Ev_f(t) + self.Im_f(t) + self.Icr_f(t) + self.Iv_f(t) + self.H_f(t) + self.D_f(t) + self.R_f(t) 


    def integrate(self,t0=0,T=None,h=0.01):
        print('The use of integrate() is now deprecated. Use solve() instead.')
        self.solve(t0=t0,T=T,h=h)

    # Scipy
    def solve(self,t0=0,T=None,h=0.01):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,mu*(self.I_ac)T] function doesn't
        #start. Or it's start if class object is initialze.

        if T is None:
            T = self.tsim

        if(not isinstance(self.S, np.ndarray)):
            S0=self.S
            Sv0=self.Sv

            if self.E:
                E0 = self.E
                E_d0 = self.E_d
                Ev0 = self.Ev
                Ev_d0 = self.Ev_d                
            else:
                E0 = self.mu*(self.I)
                E_d0 = self.mu*(self.I_d)
                
                Ev0 = self.mu*(self.Iv)
                Ev_d0 = self.mu*(self.Iv_d)


            Im0=self.Im
            Im_d0=self.Im_d
            Icr0=self.Icr
            Icr_d0=self.Icr_d
            Iv0=self.Iv
            Iv_d0=self.Iv_d

            H0 = self.H
            H_d0 = self.H_d
            
            D0 = self.D
            D_d0 = self.D_d

            R0=self.R
            R_d0=0

            Phi0=0

            self.t=np.arange(t0,T+h,h)
            

        else:
            return()
            

        initcond = np.array([S0,Sv0,E0,E_d0,Ev0,Ev_d0,Im0,Im_d0,Icr0,Icr_d0,Iv0,Iv_d0,H0,H_d0,D0,D_d0,R0,R_d0,Phi0])

        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dS(t,y[0],y[6],y[8],y[10],y[16],y[18])
            ydot[1]=self.dSv(t,y[1],y[6],y[8],y[10],y[16],y[18])

            ydot[2]=self.dE(t,y[0],y[2],y[6],y[8],y[10],y[18])
            ydot[3]=self.dE_d(t,y[0],y[3],y[6],y[8],y[10],y[18])

            ydot[4]=self.dEv(t,y[1],y[4],y[6],y[8],y[10],y[18])
            ydot[5]=self.dEv_d(t,y[1],y[5],y[6],y[8],y[10],y[18])

            ydot[6]=self.dIm(t,y[2],y[6])
            ydot[7]=self.dIm_d(t,y[2],y[7])

            ydot[8]=self.dIcr(t,y[2],y[8])
            ydot[9]=self.dIcr_d(t,y[2],y[9])

            ydot[10]=self.dIv(t,y[4],y[10])
            ydot[11]=self.dIv_d(t,y[4],y[11])                        

            ydot[12]=self.dH(t,y[8],y[10],y[12])
            ydot[13]=self.dH_d(t,y[8],y[10],y[12],y[13])

            ydot[14]=self.dD(t,y[8],y[10],y[12])
            ydot[15]=self.dD_d(t,y[8],y[10],y[12],y[15])

            ydot[16]=self.dR(t,y[6],y[10],y[12],y[16])
            ydot[17]=self.dR_d(t,y[6],y[10],y[12],y[17])

            ydot[18]=self.dPhi(t)
            return(ydot)

         
        
        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA',t_eval=list(range(t0,T)))
        
        self.sol = sol
        self.t=sol.t 
        
        self.S=sol.y[0,:]
        self.Sv=sol.y[1,:]
        self.E=sol.y[2,:]
        self.E_d=sol.y[3,:]
        self.Ev=sol.y[4,:]
        self.Ev_d=sol.y[5,:]
        self.Im=sol.y[6,:]
        self.Im_d=sol.y[7,:]
        self.Icr=sol.y[8,:]
        self.Icr_d=sol.y[9,:]
        self.Iv=sol.y[10,:]
        self.Iv_d=sol.y[11,:]
        self.H=sol.y[12,:]
        self.H_d=sol.y[13,:]
        self.D=sol.y[14,:]
        self.D_d=sol.y[15,:]
        self.R=sol.y[16,:]
        self.R_d=sol.y[17,:]
        self.Phi=sol.y[18,:]

        # Calculate accumulated variables
        self.E_ac = np.cumsum(self.E_d)
        self.Ev_ac = np.cumsum(self.Ev_d)

        self.Im_ac = np.cumsum(np.append([0],self.Im_d[:-1])) + self.Im_ac
        self.Icr_ac = np.cumsum(np.append([0],self.Icr_d[:-1])) + self.Icr_ac
        self.Iv_ac = np.cumsum(np.append([0], self.Iv_d[:-1]
        ))+ self.Iv_ac

        self.R_ac = np.cumsum(self.R_d)

        self.I = self.Im + self.Icr + self.Iv
        self.I_d = self.Im_d + self.Icr_d + self.Iv_d
        self.I_ac = self.Im_ac + self.Icr_ac + self.Iv_ac        

        self.underreport()
        self.analytics()
        self.df_build()

        return

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
        self.prevalence_susc = np.array([self.I_ac[i]/(self.S[i]+self.E[i]+self.I[i]+self.R[i]) for i in range(len(self.I_ac))])
        self.prevalence_det = np.array([self.I_ac_det[i]/(self.S[i]+self.E[i]+self.I_det[i]+self.R[i]) for i in range(len(self.I_ac))])

        self.CFR = np.array([self.pH_D(t)*self.pE_Icr(t) for t in self.t])

    def underreport(self):
        """Calculates the detected cases using the underreport factor
        """
        # ToDo: Revisar el cálculo del subreporte

        if False:
            self.Im_det = [self.Im[i]*self.pI_det(self.t(i)) for i in range(len(self.t))]
            self.Im_d_det = [self.Im_d[i]*self.pI_det(self.t(i)) for i in range(len(self.t))]
            self.Im_ac_det = np.cumsum(self.Im_d_det)

        self.Im_det = self.Im*self.pI_det(0)
        self.Im_d_det = self.Im_d*self.pI_det(0)
        self.Im_ac_det = self.Im_ac*self.pI_det(0)

        self.Icr_det = self.Icr*self.pI_det(0)
        self.Icr_d_det = self.Icr_d*self.pI_det(0)
        self.Icr_ac_det = self.Icr_ac*self.pI_det(0)

        self.Iv_det = self.Iv*self.pIv_det(0)
        self.Iv_d_det = self.Iv_d*self.pIv_det(0)
        self.Iv_ac_det = self.Iv_ac*self.pIv_det(0)


        self.I_det = self.Im_det + self.Icr_det + self.Iv_det
        self.I_d_det = self.Im_d_det + self.Icr_d_det + self.Iv_d_det
        self.I_ac_det = self.Im_ac_det + self.Icr_ac_det + self.Iv_ac_det

       
    def df_build(self):        
        self.results = pd.DataFrame({'t':self.t,'dates':self.dates})
        names = ['S','Sv','E','E_d','Ev','Ev_d','Im','Im_d','Icr','Icr_d','Iv','Iv_d','H','H_d','D','D_d','R','R_d','Phi']
        
        self.aux = pd.DataFrame(np.transpose(self.sol.y),columns=names)       

        names2 = ['E_ac','Ev_ac','Im_ac','Icr_ac','Iv_ac','R_ac','Im_det','Im_d_det','Im_ac_det','Icr_det','Icr_d_det','Icr_ac_det',
        'Iv_det','Iv_d_det','Iv_ac_det','I','I_d','I_ac','I_det','I_d_det','I_ac_det','prevalence_total','prevalence_susc','prevalence_det','CFR']
        vars2 = [self.E_ac,self.Ev_ac,self.Im_ac,self.Icr_ac,self.Iv_ac,self.R_ac,self.Im_det,self.Im_d_det,self.Im_ac_det,
        self.Icr_det,self.Icr_d_det,self.Icr_ac_det,self.Iv_det,self.Iv_d_det,self.Iv_ac_det,self.I,self.I_d,self.I_ac,
        self.I_det,self.I_d_det,self.I_ac_det,
        self.prevalence_total,self.prevalence_susc,self.prevalence_det,self.CFR]
        
        self.aux2 = pd.DataFrame(np.transpose(vars2),columns=names2)

        self.results = pd.concat([self.results,self.aux,self.aux2],axis=1)

        # Cast int
        for i in self.results.keys():
            if not i in ['dates','prevalence_total','prevalence_susc','prevalence_det','CFR']:
                self.results[i] = self.results[i].astype('int')
                #pass

        #self.results = self.results.astype({'S': int,'E': int,'E_d': int,'I': int,'I_d': int,'R': int,'R_d': int,'E_ac': int,
        #'I_ac': int,'R_ac': int,'I_det': int,'I_d_det': int,'I_ac_det': int})

        self.resume = pd.DataFrame({'peak':int(self.peak),'peak_t':self.peak_t,'peak_date':self.peak_date},index=[0])



    """

    def calculate_indicators(self):
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

        