#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from datetime import timedelta
from datetime import datetime
from numpy import linalg as LA

"""
# --------------------------------------- #
#       Simulation results processing     #
# --------------------------------------- #

Creation of local variables for later analysis
"""


class SEIRHVD_vars():

    # Creacion de variables auxiliares para los analisis
    def localvar(self):
        # Poblacion total
        self.T=[self.sims[i][0].S+self.sims[i][0].E_as+self.sims[i][0].E_sy+self.sims[i][0].I_as+self.sims[i][0].I_cr+self.sims[i][0].I_mi+self.sims[i][0].I_se\
            +self.sims[i][0].H_in+self.sims[i][0].H_out+self.sims[i][0].H_cr+self.sims[i][0].V+self.sims[i][0].D+self.sims[i][0].R+self.sims[i][0].B for i in range(self.numescenarios)]


        # Susceptibles
        self.S = [self.sims[i][0].S for i in range(self.numescenarios)]
        # Hospitalizados totales diarios
        self.H_sum=[self.sims[i][0].H_in+self.sims[i][0].H_cr+self.sims[i][0].H_out+self.sims[i][0].V for i in range(self.numescenarios)] 
        # Hospitalizados camas diarios
        self.H_bed=[self.sims[i][0].H_in+self.sims[i][0].H_cr+self.sims[i][0].H_out for i in range(self.numescenarios)] 
        # Hospitalizados ventiladores diarios
        self.H_vent=[self.sims[i][0].V for i in range(self.numescenarios)] 
        # Infectados Acumulados
        self.Iac=[self.sims[i][0].I for i in range(self.numescenarios)] 
        # Infectados activos diarios
        self.I = [self.sims[i][0].I_as+self.sims[i][0].I_cr + self.sims[i][0].I_mi + self.sims[i][0].I_se + self.sims[i][0].H_in+self.sims[i][0].H_cr+self.sims[i][0].H_out+self.sims[i][0].V for i in range(self.numescenarios)] 
        #self.I_act = [self.sims[i][0].I_mi + self.sims[i][0].I_cr + self.sims[i][0].I_se + self.sims[i][0].H_in+self.sims[i][0].H_cr+self.sims[i][0].H_out+self.sims[i][0].V for i in range(self.numescenarios)] 
        self.I_act = [self.sims[i][0].I_as+self.sims[i][0].I_mi + self.sims[i][0].I_cr + self.sims[i][0].I_se  for i in range(self.numescenarios)] 

        # Infectados asintomaticos
        self.I_as = [self.sims[i][0].I_as for i in range(self.numescenarios)]         
        # Infectados mild
        self.I_mi = [self.sims[i][0].I_mi for i in range(self.numescenarios)] 
        # Infectados severos
        self.I_se = [self.sims[i][0].I_se for i in range(self.numescenarios)] 
        # Infectados criticos
        self.I_cr = [self.sims[i][0].I_cr for i in range(self.numescenarios)] 
        # suma de infectados "sueltos"
        self.I_sum = [self.sims[i][0].I_as+self.sims[i][0].I_cr + self.sims[i][0].I_mi + self.sims[i][0].I_se for i in range(self.numescenarios)] 

        # Infectados nuevos diarios
        # Infectados asintomaticos
        self.I_as_d = [self.sims[i][0].I_as_d for i in range(self.numescenarios)]         
        # Infectados mild
        self.I_mi_d = [self.sims[i][0].I_mi_d for i in range(self.numescenarios)] 
        # Infectados severos
        self.I_se_d = [self.sims[i][0].I_se_d for i in range(self.numescenarios)] 
        # Infectados criticos
        self.I_cr_d = [self.sims[i][0].I_cr_d for i in range(self.numescenarios)] 


        # Infectados acumulados
        # Infectados asintomaticos
        self.I_as_ac = [self.sims[i][0].I_as_ac for i in range(self.numescenarios)]         
        # Infectados mild
        self.I_mi_ac = [self.sims[i][0].I_mi_ac for i in range(self.numescenarios)] 
        # Infectados severos
        self.I_se_ac = [self.sims[i][0].I_se_ac for i in range(self.numescenarios)] 
        # Infectados criticos
        self.I_cr_ac = [self.sims[i][0].I_cr_ac for i in range(self.numescenarios)] 


        # Expuestos totales diarios
        self.E = [self.sims[i][0].E_as+self.sims[i][0].E_sy for i in range(self.numescenarios)]
        self.E_as = [self.sims[i][0].E_as for i in range(self.numescenarios)]  
        self.E_sy = [self.sims[i][0].E_sy for i in range(self.numescenarios)]  

        # Fallecidos::

        # Enterrados/Muertos acumulados
        self.B = [self.sims[i][0].B for i in range(self.numescenarios)] 
        # Muertos diarios
        self.D = [self.sims[i][0].D for i in range(self.numescenarios)] 
        # Infectados Criticos Fallecidos
        self.I_crD = [self.sims[i][0].I_crD for i in range(self.numescenarios)] 
        self.I_crD_d = [self.sims[i][0].I_crD_d for i in range(self.numescenarios)]
        # Infectados Serios Fallecidos
        self.I_seD = [self.sims[i][0].I_seD for i in range(self.numescenarios)]
        self.I_seD_d = [self.sims[i][0].I_seD_d for i in range(self.numescenarios)]
        # Hospitalizados Criticos Fallecidos
        self.H_crD = [self.sims[i][0].H_crD for i in range(self.numescenarios)]
        self.H_crD_d = [self.sims[i][0].H_crD_d for i in range(self.numescenarios)]
        # Ventilados Fallecidos
        self.VD = [self.sims[i][0].VD for i in range(self.numescenarios)]
        self.VD_d = [self.sims[i][0].VD_d for i in range(self.numescenarios)]

        # Recuperados
        self.R = [self.sims[i][0].R for i in range(self.numescenarios)] 
        # Ventiladores diarios
        self.V = [self.sims[i][0].V for i in range(self.numescenarios)] 

        # Variables temporales
        self.t = [self.sims[i][0].t for i in range(self.numescenarios)] 
        self.dt = [np.diff(self.t[i]) for i in range(self.numescenarios)] 
        
        
        # CAMAS
        self.H_crin=[self.sims[i][0].H_cr for i in range(self.numescenarios)] 
        self.H_in=[self.sims[i][0].H_in for i in range(self.numescenarios)] 
        self.H_out=[self.sims[i][0].H_out for i in range(self.numescenarios)] 
        self.H_sum=[self.sims[i][0].H_in+self.sims[i][0].H_cr+self.sims[i][0].H_out for i in range(self.numescenarios)]
        self.H_tot=[self.sims[i][0].H_in+self.sims[i][0].H_cr+self.sims[i][0].H_out+self.sims[i][0].V  for i in range(self.numescenarios)]

        self.CH = [self.sims[i][0].CH for i in range(self.numescenarios)]
        self.CV = [self.sims[i][0].CV for i in range(self.numescenarios)]
        self.ACH = [self.sims[i][0].ACH for i in range(self.numescenarios)]
        self.ACV = [self.sims[i][0].ACV for i in range(self.numescenarios)]
        
        #Cálculo de la fecha del Peak  
        self.peakindex = [np.where(self.I[i]==max(self.I[i]))[0][0] for i in range((self.numescenarios))]
        self.peak = [max(self.I[i]) for i in range((self.numescenarios))]
        self.peak_t = [self.t[i][self.peakindex[i]] for i in range((self.numescenarios))]
        self.peak_date = [self.initdate+timedelta(days=round(self.peak_t[i])) for i in range((self.numescenarios))]

        #Cálculo de la fecha de Saturacion_
        self.ventsat = 0
        self.hsta = 0

        #proporcion de la poblacion que entra en la dinamica de infeccion
        #self.population = self.sims[0][0].population
        self.infectedsusc = [100*((self.S[i][0] - self.S[i][-1])/self.S[i][0]) for i in range(self.numescenarios)] 
        self.infectedpop = [100*((self.S[i][0] - self.S[i][-1]))/self.population for i in range(self.numescenarios)]


        # Indicadores:
        # Should Have been hospitalized
        #self.SHFR = [self.B[i]/(self.I_se_ac[i]+self.I_cr_ac[i]) for i in range(self.numescenarios)]
        self.totD =  [self.B[i][-1] for i in range(self.numescenarios)]
        
        #self.SHFR = [self.totD[i]/(self.I_se_ac[i][-1]+self.I_cr_ac[i][-1]) if self.totD[i]>15 else 0.15 for i in range(self.numescenarios)] 
        self.SHFR = [self.totD[i]/(self.I_se_ac[i][-1]+self.I_cr_ac[i][-1]) for i in range(self.numescenarios)] 
        self.SHFR_d = [self.B[i]/(self.I_se_ac[i][-1]+self.I_cr_ac[i][-1]) for i in range(self.numescenarios)]

        # ----------------- #
        #    QA Variables   #
        # ----------------- #
        # Infected accumulated checksum
        self.Iac_checksum = [[self.Iac[i][j] - self.I_as_ac[i][j] - self.I_mi_ac[i][j] - self.I_se_ac[i][j] - self.I_cr_ac[i][j] - self.I_act0 for j in range(len(self.t[i]))] for i in range(self.numescenarios)]
        # Accumulated Infected proportion over time
        self.I_as_ac_prop = [[self.I_as_ac[i][j]/self.Iac[i][j] for j in range(len(self.t[i]))] for i in range(self.numescenarios)]
        self.I_mi_ac_prop = [[self.I_mi_ac[i][j]/self.Iac[i][j] for j in range(len(self.t[i]))] for i in range(self.numescenarios)]
        self.I_se_ac_prop = [[self.I_se_ac[i][j]/self.Iac[i][j] for j in range(len(self.t[i]))] for i in range(self.numescenarios)]
        self.I_cr_ac_prop = [[self.I_cr_ac[i][j]/self.Iac[i][j] for j in range(len(self.t[i]))] for i in range(self.numescenarios)]
 
        # Accumulated Infected proportion sum - must be 1
        self.Iac_prop = [self.I_as_ac_prop[i][-1] +self.I_mi_ac_prop[i][-1] +self.I_se_ac_prop[i][-1] +self.I_cr_ac_prop[i][-1] +self.I_act0/self.Iac[i][-1] for i in range(self.numescenarios)]

        # Deaths
        self.D_checksum = [max(self.VD_d[i]+self.H_crD_d[i]+self.I_seD_d[i]+self.I_crD_d[i] -self.D[i]) for i in range(self.numescenarios)]
        #self.B_checksum = [max(self.VD[i]+self.H_crD[i]+self.I_seD[i]+self.I_crD[i]-self.B[i]) for i in range(self.numescenarios)]

        # Population:
        self.population_checksum = [max(self.S[i]+self.E_as[i]+self.E_sy[i]+self.I_as[i]+self.I_mi[i]+self.I_se[i]+self.I_cr[i]
            +self.H_crin[i]+self.H_in[i]+self.H_out[i]+self.V[i]+self.R[i]+self.D[i]-self.population) for i in range(self.numescenarios)]

        # -------------- #
        #     Errores    #
        # -------------- #
        if self.realdata:
            # Camas
            idx = [np.searchsorted(self.t[i],self.sochimi_tr) for i in range(self.numescenarios)]
            self.err_bed = [LA.norm(self.Hr-self.H_sum[i][idx[i]])/LA.norm(self.Hr) for i in range(self.numescenarios)]
            self.err_vent = [LA.norm(self.Vr-self.V[i][idx[i]])/LA.norm(self.Vr) for i in range(self.numescenarios)]  
            
            # Infecatos Activos
            idx = [np.searchsorted(self.t[i],self.tr) for i in range(self.numescenarios)]
            self.err_Iactives = [LA.norm(self.Ir-self.I[i][idx[i]])/LA.norm(self.Ir) for i in range(self.numescenarios)]    
            
            # Infectados acumulados
            #idx = [np.searchsorted(t[i],tr) for i in range(self.numescenarios)]
            #err_Iactives = [LA.norm(Ir-I[i][idx[i]])/LA.norm(Ir) for i in range(self.numescenarios)]    
            
            # Fallecidos
            idx = [np.searchsorted(self.t[i],self.Br_tr) for i in range(self.numescenarios)]
            self.err_dead = [LA.norm(self.Br-self.B[i][idx[i]])/LA.norm(self.Br) for i in range(self.numescenarios)]

        self.quarantines = []
        for i in range(self.numescenarios):
            self.quarantines.append([self.sims[i][0].alpha(t) for t in self.t[i]])        

        # -------------------- #
        #   Fecha de Colapso   #
        # -------------------- #
        self.H_colapsedate = []
        self.H_colapseday = []
        self.V_colapsedate = []
        self.V_colapseday = []
        for i in range(self.numescenarios):
            try:
                self.H_colapseday.append(np.where(self.CH[i]>0)[0][0])
                self.H_colapsedate.append((self.initdate+timedelta(days=self.H_colapseday[i])).strftime('%Y/%m/%d'))
            except:
                self.H_colapseday.append(self.tsim)
                self.H_colapsedate.append("None")
            try:
                self.V_colapseday.append(np.where(self.CV[i]>0)[0][0])
                self.V_colapsedate.append((self.initdate+timedelta(days=self.V_colapseday[i])).strftime('%Y/%m/%d'))
            except:
                self.V_colapseday.append(self.tsim)
                self.V_colapsedate.append("None")

    """
        # ---------------------------------- #
        #          Estudio Resultados        #
        # ---------------------------------- # 
    """
    # ---------------------------------- #
    #        Resumen de resultados       #
    # ---------------------------------- #
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
