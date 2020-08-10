#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
   
   SEIR Plot Functions

"""
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA


class SEIR_plots():
    # -------------------------- #
    #        Plot function       #
    # -------------------------- #
    def plot(self,title = '',xlabel='',ylabel='',legend=True):
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        if legend:
            plt.legend(loc=0)
        plt.show()

    # ---------------------------------------- #
    #                 Datos                    #
    # ---------------------------------------- #
    def plotdatosactivos(self,enddate =  datetime(2020,7,30),days=0, reales= False,ylim = 0,norm=1,scalefactor = False,legend=True):
        if not self.realdata:
            return('No real data')        
        # Reales
        if reales:
            plt.scatter(self.tr,self.Ir,label='Infectados Activos reales')

        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
          
        self.plot(title = 'Activos Reales',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)

    def plotdatosacumulados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True):
        if not self.realdata:
            return('No real data')        
        # Reales
        if reales:
            plt.scatter(self.I_ac_r_tr,self.I_ac_r,label='Infectados Acumulados reales')

        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
          
        self.plot(title = 'Infectados Acumulados Reales - EPI',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)

    def plotdatosdiarios(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True):
        if not self.realdata:
            return('No real data')        
        # Reales
        if reales:
            plt.scatter(self.I_d_r_tr,self.I_d_r,label='Infectados diarios reales')

        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
          
        self.plot(title = 'Infectados Diarios Reales - EPI',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)



    # -------------------------------------------------------- #
    #                       Infectados                         #
    # -------------------------------------------------------- #

    # ------------------------------ #
    #       Infectados Activos       #
    # ------------------------------ #
    def plotActiveInfected(self,enddate =  datetime(2020,7,30),days=-1, reales= False,ylim = 0,norm=1,scalefactor = False,legend=True, minciencia = True,showparams=False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days
        if days < 0:
            days = self.tsim
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]

        Isf = norm    
        if scalefactor:
            Isf = self.ScaleFactor*norm


        # ----------- #
        #     Plot    #
        # ----------- #
        
        # Error
        if self.realdata:                        
            if minciencia:
                for i in range(self.numescenarios):
                    plt.plot([], [], ' ', label='err: '+str(round(100*self.err_Iactives_minciencia[i],2))+'%'+' Mov = '+str(self.inputarray[i][2]))    
        
            else:
                for i in range(self.numescenarios):
                    plt.plot([], [], ' ', label='err: '+str(round(100*self.err_Iactives[i],2))+'%'+' Mov = '+str(self.inputarray[i][2]))    

        # Parametros:
        if showparams:
            plt.plot([], [], ' ', label='beta: '+str(self.beta))
            plt.plot([], [], ' ', label='mu: '+str(self.mu))
            plt.plot([], [], ' ', label='k: '+str(self.k))
            #plt.plot([], [], ' ', label='ScaleFactor: '+str(self.ScaleFactor))
            #plt.plot([], [], ' ', label='seroprev: '+str(self.SeroPrevFactor))
            if scalefactor:
                for i in range(self.numescenarios):
                    plt.plot([], [], ' ', label='SeroPrev Norm: '+str(round(self.infectedpop_norm[i],2))+'%'+' Mov = '+str(self.inputarray[i][2]))    
            else:
                for i in range(self.numescenarios):
                    plt.plot([], [], ' ', label='SeroPrev: '+str(round(self.infectedpop[i],2))+'%'+' Mov = '+str(self.inputarray[i][2]))    


        # Reales
        if self.realdata:
            if reales:
                if minciencia:
                    plt.scatter(self.I_minciencia_r_tr,self.I_minciencia_r,label='Infectados Activos Minciencia')
                else:
                    plt.scatter(self.tr,self.Ir,label='Infectados Activos ')                    


        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dotted',color = 'grey')

        # Infectados
        #linestyle = ['dashed','solid','dotted','dashed','solid','dotted','solid','dashed','dotted','dotted']
        #linestyle = ['dashed','solid','dashed','solid','dotted']
        colors = ['red','blue','green','purple','black','lime','cyan','m','indigo','orange','orangered','wheat','salmon']
        #colors = ['lime','lime','purple','purple','black']
        qt = ['TQ' if self.inputarray[i][-1]==0 else 'DQ'+str(int(self.inputarray[i][3])) for i in range(len(self.inputarray))]
        for i in range(self.numescenarios):        
            plt.plot(self.t[i],self.I_act[i]/Isf,label='Mob = '+str(self.inputarray[i][2])+' '+qt[i])#,color = colors[i],linestyle='solid',linewidth=2)

        if days >0:
            plt.xlim(0,days)
        if ylim >0:
            plt.ylim(0,ylim)            
        self.plot(title = 'Active infected',xlabel='Day',legend=legend)


  
    # -------------------------------- #
    #       Infectados Acumulados      #
    # -------------------------------- #
    def plotAccumulatedInfected(self,enddate =  datetime(2020,7,30),days=-1, reales= False,ylim = 0,norm=1,scalefactor = False,showparams = False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days
        if days < 0:
            days = self.tsim
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        
        Isf = norm    
        if scalefactor:
            Isf = self.ScaleFactor*norm

        # ----------- #
        #     Plot    #
        # ----------- #
        if showparams:
            plt.plot([], [], ' ', label='beta: '+str(self.beta))
            plt.plot([], [], ' ', label='mu: '+str(self.mu))
            plt.plot([], [], ' ', label='k: '+str(self.k))        
        
        # Error
        #if showerror:
        #    for i in range(self.numescenarios):
        #        plt.plot([], [], ' ', label='err: '+str(round(100*self.err[i],2))+'%')    

        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dotted',color = 'grey')        

        # Reales
        if self.realdata:
            if reales:
                plt.scatter(self.I_ac_r_tr,self.I_ac_r,label='Infectados Acumulados reales')

        # Infectados
        qt = ['TQ' if self.inputarray[i][-1]==0 else 'DQ'+str(int(self.inputarray[i][3])) for i in range(len(self.inputarray))]
        for i in range(self.numescenarios):        
            plt.plot(self.t[i],self.Iac[i]/Isf,label='Mob = '+str(self.inputarray[i][2])+' '+qt[i])

        if days >0:
            plt.xlim(0,days)
        if ylim >0:
            plt.ylim(0,ylim)            
        self.plot(title = 'Accumulated Infected',xlabel='Day')  


    # ------------------------------ #
    #       Infectados Diarios       #
    # ------------------------------ #    
    def plotDailyInfected(self,enddate =  datetime(2020,7,30),days=-1, reales= False,ylim = 0,norm=1,scalefactor = False, showparams=False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days
        if days < 0:
            days = self.tsim
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)] 

        idx = [np.searchsorted(self.t[i],range(days)) for i in range(self.numescenarios)]

        # Reales
        if self.realdata:
            if reales:
                plt.scatter(self.I_d_r_tr,self.I_d_r,label='Infectados diarios reales')

        if showparams:
            plt.plot([], [], ' ', label='beta: '+str(self.beta))
            plt.plot([], [], ' ', label='mu: '+str(self.mu))
            plt.plot([], [], ' ', label='k: '+str(self.k))

        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')       

        Idiarios = [self.I_as_d[i]+self.I_mi_d[i]+self.I_se_d[i]+self.I_cr_d[i] for i in range(self.numescenarios)]

       
        qt = ['TQ' if self.inputarray[i][-1]==0 else 'DQ'+str(int(self.inputarray[i][3])) for i in range(len(self.inputarray))]
        for i in range(self.numescenarios):        
            plt.plot(self.t[i][:(endD[i])],Idiarios[i][:endD[i]],label='Mob = '+str(self.inputarray[i][2])+' '+qt[i])

        if days >0:
            plt.xlim(0,days)
        self.plot(title = 'Infectados Diarios',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))          
  


    # -------------------- #
    #     Susceptible      #
    # -------------------- #
    def plotSusceptible(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True,showparams=False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        
        if norm <1:
            norm = self.ScaleFactor

        if showparams:
            plt.plot([], [], ' ', label='beta: '+str(self.beta))
            plt.plot([], [], ' ', label='mu: '+str(self.mu))
            plt.plot([], [], ' ', label='k: '+str(self.k))            
        #Isf = 1    
        #if scalefactor:
        #    Isf = ScaleFactor

        # ----------- #
        #     Plot    #
        # ----------- #
        # Parametros 
        #plt.plot([], [], ' ', label='beta: '+str(self.beta))
        #plt.plot([], [], ' ', label='mu: '+str(self.mu))
        #plt.plot([], [], ' ', label='factor de escala: '+str(self.ScaleFactor))

        # Fecha de Peak
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        #

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        qt = ['TQ' if self.inputarray[i][-1]==0 else 'DQ'+str(int(self.inputarray[i][3])) for i in range(len(self.inputarray))]
        for i in range(self.numescenarios):
            plt.plot(self.t[i][:endD[i]],self.S[i][:endD[i]],label='Mob = '+str(self.inputarray[i][2])+' '+qt[i])#,color = 'blue',linestyle=linestyle[i])
            #plt.plot(self.t[i][:endD[i]],self.E_sy[i][:endD[i]],label='Expuestos sintomáticos Mov = '+str(self.inputarray[i][2]),color = 'red',linestyle=linestyle[i])
            
        plt.xlim(0,days)   
        if ylim >0:
            plt.ylim(0,ylim)

        self.plot(title = 'Susceptible',xlabel='Days',legend=legend)
        

    # ------------------ #
    #     Expuestos      #
    # ------------------ #
    def plotExposed(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True,showparams=False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        
        if norm <1:
            norm = self.ScaleFactor

        if showparams:
            plt.plot([], [], ' ', label='beta: '+str(self.beta))
            plt.plot([], [], ' ', label='mu: '+str(self.mu))
            plt.plot([], [], ' ', label='k: '+str(self.k))            
        #Isf = 1    
        #if scalefactor:
        #    Isf = ScaleFactor

        # ----------- #
        #     Plot    #
        # ----------- #
        # Parametros 
        #plt.plot([], [], ' ', label='beta: '+str(self.beta))
        #plt.plot([], [], ' ', label='mu: '+str(self.mu))
        #plt.plot([], [], ' ', label='factor de escala: '+str(self.ScaleFactor))

        # Fecha de Peak
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        #

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        qt = ['TQ' if self.inputarray[i][-1]==0 else 'DQ'+str(int(self.inputarray[i][3])) for i in range(len(self.inputarray))]
        for i in range(self.numescenarios):
            plt.plot(self.t[i][:endD[i]],self.E[i][:endD[i]],label='Mob = '+str(self.inputarray[i][2])+' '+qt[i])#,color = 'blue',linestyle=linestyle[i])
            #plt.plot(self.t[i][:endD[i]],self.E_sy[i][:endD[i]],label='Expuestos sintomáticos Mov = '+str(self.inputarray[i][2]),color = 'red',linestyle=linestyle[i])
            
        plt.xlim(0,days)   
        if ylim >0:
            plt.ylim(0,ylim)

        self.plot(title = 'Exposed',xlabel='Days',legend=legend)
        

    # ------------------ #
    #     Recovered      #
    # ------------------ #
    def plotRecovered(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True,showparams=False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        
        if norm <1:
            norm = self.ScaleFactor

        if showparams:
            plt.plot([], [], ' ', label='beta: '+str(self.beta))
            plt.plot([], [], ' ', label='mu: '+str(self.mu))
            plt.plot([], [], ' ', label='k: '+str(self.k))            
        #Isf = 1    
        #if scalefactor:
        #    Isf = ScaleFactor

        # ----------- #
        #     Plot    #
        # ----------- #
        # Parametros 
        #plt.plot([], [], ' ', label='beta: '+str(self.beta))
        #plt.plot([], [], ' ', label='mu: '+str(self.mu))
        #plt.plot([], [], ' ', label='factor de escala: '+str(self.ScaleFactor))

        # Fecha de Peak
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        #

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        qt = ['TQ' if self.inputarray[i][-1]==0 else 'DQ'+str(int(self.inputarray[i][3])) for i in range(len(self.inputarray))]
        for i in range(self.numescenarios):
            plt.plot(self.t[i][:endD[i]],self.R[i][:endD[i]],label='Mob = '+str(self.inputarray[i][2])+' '+qt[i])#,color = 'blue',linestyle=linestyle[i])
            #plt.plot(self.t[i][:endD[i]],self.E_sy[i][:endD[i]],label='Expuestos sintomáticos Mov = '+str(self.inputarray[i][2]),color = 'red',linestyle=linestyle[i])
            
        plt.xlim(0,days)   
        if ylim >0:
            plt.ylim(0,ylim)

        self.plot(title = 'Recovered',xlabel='Days',legend=legend)
        

    # -------------------- #
    #     Curvas SEIR      #
    # -------------------- #
    def plotseir(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,seir = [1,1,1,1],showparams=False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        

        if norm <1:
            norm = self.ScaleFactor
        #Isf = 1    
        #if scalefactor:
        #    Isf = ScaleFactor

        # ----------- #
        #     Plot    #
        # ----------- #
        # Parametros 
        if showparams:
            plt.plot([], [], ' ', label='beta: '+str(self.beta))
            plt.plot([], [], ' ', label='mu: '+str(self.mu))
            plt.plot([], [], ' ', label='k: '+str(self.k))         

        # Fecha de Peak
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        

        linestyle = ['solid','dashed','dotted','-.']
        colors = ['red','blue','green','purple','black','lime','cyan','m','indigo','orange','orangered','wheat','salmon']
        qt = ['TQ' if self.inputarray[i][-1]==0 else 'DQ'+str(int(self.inputarray[i][3])) for i in range(len(self.inputarray))]        
        for i in range(self.numescenarios):        
            plt.plot(self.t[i],self.S[i],label='S Mob = '+str(self.inputarray[i][2])+' '+qt[i],linestyle=linestyle[0],color = colors[i])
            plt.plot(self.t[i],self.E[i],label='E Mob = '+str(self.inputarray[i][2])+' '+qt[i],linestyle=linestyle[2],color = colors[i])
            plt.plot(self.t[i],self.I[i],label='I Mob = '+str(self.inputarray[i][2])+' '+qt[i],linestyle=linestyle[1],color = colors[i])        
            plt.plot(self.t[i],self.R[i],label='R Mob = '+str(self.inputarray[i][2])+' '+qt[i],linestyle=linestyle[3],color = colors[i])
            #plt.plot(self.t[i],D[i],label='Muertos diarios Mov = '+str(inputarray[i][2]),linestyle=linestyle[i])
            #plt.plot(self.t[i],self.B[i],label='Enterrados Mov = '+str(self.inputarray[i][2]),linestyle=linestyle[i])
            
        plt.xlim(0,days)   
        if ylim >0:
            plt.ylim(0,ylim)

        self.plot(title = 'SEIR Model',xlabel='Days')
        

    # ------------------------- #
    #       Plot Cuarentenas    #
    # ------------------------- #
    def plotQuarantines(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days
        if days < 0:
            days = self.tsim
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        
        # ----------- #
        #     Plot    #
        # ----------- #
        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dotted',color = 'grey')


        # cuarentenas
        linestyle = ['dashed','dotted','-.','*',':']
        qt = ['TQ' if self.inputarray[i][-1]==0 else 'DQ'+str(int(self.inputarray[i][3])) for i in range(len(self.inputarray))]
        for i in range(self.numescenarios):        
            plt.plot(self.t[i],self.quarantines[i],label='Mob = '+str(self.inputarray[i][2])+' '+qt[i])

        if days >0:
            plt.xlim(0,days)
        self.plot(title = 'Quarantines',xlabel='Days', ylabel='Mobility')  
