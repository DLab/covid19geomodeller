#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
   
   SEIRHDV Plot Functions

"""
import matplotlib.pyplot as plt
from numpy import linalg as LA
import numpy as np
from datetime import datetime
from datetime import timedelta

class SEIRHVD_plots():
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
    def plotdatossochimi(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        if not self.realdata:
            return('No real data')
        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')

        #ploteo de datos reales
        plt.scatter(self.sochimi_tr,self.Hr,label='Camas Ocupadas reales')
        plt.scatter(self.sochimi_tr,self.Vr,label='Ventiladores Ocupados reales')
        plt.scatter(self.sochimi_tr,self.Hr_tot,label='Capacidad de Camas')
        plt.scatter(self.sochimi_tr,self.Vr_tot,label='Capacidad de Ventiladores')
        
        self.plot(title = 'Datos Sochimi',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))        

    def plotdatosventiladores(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        if not self.realdata:
            return('No real data')        
        # -------- #
        #   Time   #
        # -------- #
        days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]


        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')

        #ploteo de datos reales
        plt.scatter(self.sochimi_tr,self.Vr,label='Ventiladores Ocupados reales')
        plt.scatter(self.sochimi_tr,self.Vr_tot,label='Capacidad de Ventiladores')
        
        self.plot(title = 'Datos Sochimi - Ventiladores',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))   

    def plotdatoscamas(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        if not self.realdata:
            return('No real data')        
        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')

        #ploteo de datos reales
        plt.scatter(self.sochimi_tr,self.Hr,label='Camas Ocupadas reales')
        plt.scatter(self.sochimi_tr,self.Hr_tot,label='Capacidad de Camas')
        
        self.plot(title = 'Datos Sochimi - Camas UCI/UTI',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))        

    def plotdatosactivos(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True):
        if not self.realdata:
            return('No real data')        
        # Reales
        if reales:
            plt.scatter(self.tr,self.Ir,label='Infectados Activos reales')

        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
          
        self.plot(title = 'Activos Reales',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)




    def plotdatosfallecidosacumulados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True):
        if not self.realdata:
            return('No real data')        
        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')

        # Datos reales
        if reales:
            plt.scatter(self.Br_tr,self.Br,label='Fallecidos reales')
            #plt.scatter(self.ED_tr,self.ED_RM_ac,label='Fallecidos excesivos proyectados')

        self.plot(title = 'Fallecidos',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)

    def plotdatosfallecidosexcesivos(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        if not self.realdata:
            return('No real data')        
        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')

        # Datos reales                 
        plt.scatter(self.ED_tr,self.ED_RM_ac,label='Fallecidos excesivos proyectados')
        self.plot(title = 'Fallecidos Excesivos',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)
    # -------------------------------------------------------- #
    #                  Uso Hospitalario                        #
    # -------------------------------------------------------- #

    # -------------------------------------- #
    #       Hospitalizados desagregados      #
    # -------------------------------------- #
    # Hospitalizados desagregados
    def plothospitalizados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        # -------- #
        #   Time   #
        # -------- #
        days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]

        for i in range(self.numescenarios):
            plt.plot(self.t[i][:endD[i]],self.H_in[i][:endD[i]],label='Hin',linestyle = 'solid')
            plt.plot(self.t[i][:endD[i]],self.H_out[i][:endD[i]],label='Hout',linestyle = 'solid')
            plt.plot(self.t[i][:endD[i]],self.H_crin[i][:endD[i]],label='Hcr_in',linestyle = 'solid')

        self.plot(title = 'Hospitalizados',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))



    # ------------------ #
    #     Ventiladores   #
    # ------------------ #
    def plotventiladores(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        # -------- #
        #   Time   #
        # -------- # 
        days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]


        # Inicio cuarentena general

        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dotted',color = 'grey')

        # Ploteo datos reales
        if self.realdata:
            if reales:            
                plt.scatter(self.sochimi_tr,self.Vr,label='Ventiladores Ocupados reales')
                plt.scatter(self.sochimi_tr,self.Vr_tot,label='Capacidad de Ventiladores')

        
        # Error y parámetros
        if self.realdata:
            for i in range(self.numescenarios):        
                plt.plot([], [], ' ', label='err_vent: '+str(round(100*self.err_vent[i],2))+'%')

        plt.plot([], [], ' ', label='beta: '+str(self.beta))
        plt.plot([], [], ' ', label='mu: '+str(self.mu))
        plt.plot([], [], ' ', label='factor B-Y: '+str(self.ScaleFactor))

        # Fecha de peaks
        for i in range(self.numescenarios):        
            plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        

        # funcion de ventiladores totales
        Vtot = [self.sims[0][0].Vtot(i) for i in self.t[0][:endD[0]]]    
        plt.plot(self.t[0][:endD[0]],Vtot,color='lime')

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        for i in range(self.numescenarios):            
            plt.plot(self.t[i][:endD[i]],self.V[i][:endD[i]],label='VMI Utilizados mov='+str(self.inputarray[i][2]),color = 'blue' ,linestyle = linestyle[i])
    

        plt.xlim(0,days)
        self.plot(title = 'Ventiladores',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))



    # ------------ #
    #     Camas    #
    # ------------ #
    def plotcamas(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        # -------- #
        #   Time   #
        # -------- #
        days = (enddate-self.initdate).days    
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        #idx = [np.searchsorted(self.t[i],self.sochimi_tr) for i in range(self.numescenarios)]

        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dotted',color = 'grey')

        # Ploteo datos reales
        if self.realdata:
            if reales:         
                plt.scatter(self.sochimi_tr,self.Hr,label='Camas Ocupadas reales')
                plt.scatter(self.sochimi_tr,self.Hr_tot,label='Capacidad de Camas')


        # Display de Parametros y errores
        if self.realdata:
            for i in range(self.numescenarios):
                plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'err_bed: '+str(round(100*self.err_bed[i],2))+'%')            
        for i in range(self.numescenarios):            
            plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        plt.plot([], [], ' ', label='beta: '+str(self.beta))
        plt.plot([], [], ' ', label='mu: '+str(self.mu))
        plt.plot([], [], ' ', label='fScale: '+str(self.ScaleFactor))
            
        
        # funcion de camas totales
        Htot = [self.sims[0][0].Htot(i) for i in self.t[0][:endD[0]]]
        plt.plot(self.t[0][:endD[0]],Htot,color='lime')
        
        for i in range(self.numescenarios):
            plt.plot(self.t[i][:endD[i]],self.H_bed[i][:endD[i]],label='Camas utilizadas mov='+str(self.inputarray[i][2]),color = 'red' ,linestyle = 'dashed')
        
        #plt.plot(self.t[i][:endD[i]],self.H_bed[i][:endD[i]],label='Camas utilizadas mov='+str(inputarray[i][2]),color = 'red' ,linestyle = 'solid')
        
        plt.xlim(0,days)
        self.plot(title = 'Camas',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))



    # -------------------------- #
    #     Camas y Ventiladores   #
    # -------------------------  #
    def plotcamasyventiladores(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days
        if days < 0:
            days = self.tsim
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]


        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dotted',color = 'grey')

        #ploteo de datos reales
        if self.realdata:
            if reales:         
                plt.scatter(self.sochimi_tr,self.Hr,label='Camas Ocupadas reales')
                plt.scatter(self.sochimi_tr,self.Vr,label='Ventiladores Ocupados reales')
                plt.scatter(self.sochimi_tr,self.Hr_tot,label='Capacidad de Camas')
                plt.scatter(self.sochimi_tr,self.Vr_tot,label='Capacidad de Ventiladores')

        if self.realdata:
            for i in range(self.numescenarios):
                plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+' err_bed: '+str(round(100*self.err_bed[i],2))+'%')

        for i in range(self.numescenarios):            
            plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        plt.plot([], [], ' ', label='beta: '+str(self.beta))
        plt.plot([], [], ' ', label='mu: '+str(self.mu))
        plt.plot([], [], ' ', label='factor de escala: '+str(self.ScaleFactor))
    

        # Camas y ventiladores totales
        Htot = [self.sims[0][0].Htot(i) for i in self.t[0][:endD[0]]]
        Vtot = [self.sims[0][0].Vtot(i) for i in self.t[0][:endD[0]]]
        plt.plot(self.t[0][:endD[0]],Htot,color='lime')
        plt.plot(self.t[0][:endD[0]],Vtot,color='lime')

        
        for i in range(self.numescenarios): 
            plt.plot(self.t[i][:endD[i]],self.H_bed[i][:endD[i]],label='Camas utilizadas mov='+str(self.inputarray[i][2]),color = 'red' ,linestyle = 'dashed')
            plt.plot(self.t[i][:endD[i]],self.V[i][:endD[i]],label='VMI Utilizados mov='+str(self.inputarray[i][2]),color = 'blue' ,linestyle = 'dashed')
        #plt.plot(self.t[i][:endD[i]],H_crin[i][:endD[i]],label='Camas críticas mov='+str(inputarray[i][2]),color = 'black' ,linestyle = 'dashed')
        
        plt.xlim(0,days)
        self.plot(title = 'Camas',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))


    # ---------- #
    #    Hrate   #
    # ---------- #
    def plothrate(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        # -------- #
        #   Time   #
        # -------- #
        days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]

        
        Hrate = [self.H_in[i]/self.H_out[i] for i in range(self.numescenarios)]

        xlabel = 'Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d')
        fig, axs = plt.subplots(3)
        #fig.suptitle(title)
        
        
        colors = ['red','blue','green','lime']
        
        
        for i in range(self.numescenarios):
            axs[0].plot(self.t[i][:endD[i]],Hrate[i][:endD[i]], label='Mov='+str(self.inputarray[i][2]),linestyle = 'solid',color = colors[i])
            axs[0].legend()
        
        
        for i in range(self.numescenarios):
            axs[1].plot(self.t[i][:endD[i]],(self.H_in[i]/self.H_sum[i])[:endD[i]],label='Hin '+'Mov='+str(self.inputarray[i][2]),linestyle = 'solid',color = colors[i])        
            axs[1].legend()        
    
        for i in range(self.numescenarios):        
            axs[2].plot(self.t[i][:endD[i]],(self.H_out[i]/self.H_sum[i])[:endD[i]],label='Hout '+'Mov='+str(self.inputarray[i][2]),linestyle = 'solid',color = colors[i])
            axs[2].legend()

        axs[0].set_title('Rate Hin/Hout')
        axs[1].set_title('Rate Hin/ Hsum')
        axs[1].set_title('Rate Hout/Hsum')
        for ax in axs.flat:
            ax.label_outer()    
        plt.xlabel(xlabel)
        plt.xlim=days                
        plt.show()
    

    # --------------------------- #
    #      Camas requeridas       #
    # --------------------------- #
    def plotcamasrequeridas(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        # ----------- #
        #     Time    #
        # ----------- #    
        if days ==0:
            days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        
        # ----------- #
        #     Plot    #
        # ----------- #
        # Fechas de colapso
        i=1
        plt.plot([], [], ' ', label='Fecha colapso Camas: '+str(round(self.t[i][self.H_colapsedate[i]])))
        plt.plot([], [], ' ', label='Fecha colapso Ventiladores: '+str(round(self.t[i][self.V_colapsedate[i]])))

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        for i in range(self.numescenarios):
            plt.plot(self.t[i][:endD[i]],self.CH[i][:endD[i]],label='Intermedio/Intensivo Mov = '+str(self.inputarray[i][1]),color = 'red' ,linestyle = linestyle[i])
            plt.plot(self.t[i][:endD[i]],self.CV[i][:endD[i]],label='VMI Mov = '+str(self.inputarray[i][1]),color = 'blue' ,linestyle = linestyle[i])
        
        self.plot(title='Camas Requeridas',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))


    # ------------------------------------ #
    #      Necesidad total de Camas        #
    # ------------------------------------ #
    def plotnecesidadtotcamas(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
        # ----------- #
        #     Time    #
        # ----------- #    
        if days ==0:
            days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]

        # ----------- #
        #     Plot    #
        # ----------- #
        
        # Fechas de colapso
        i=1
        plt.plot([], [], ' ', label='Fecha colapso Camas: '+str(round(self.t[i][self.H_colapsedate[i]])))
        plt.plot([], [], ' ', label='Fecha colapso Ventiladores: '+str(round(self.t[i][self.V_colapsedate[i]])))

        # Datos reales
        if self.realdata:
            if reales:         
                plt.scatter(self.sochimi_tr,self.Hr,label='Camas Ocupadas reales')
                plt.scatter(self.sochimi_tr,self.Vr,label='Ventiladores Ocupados reales')
                plt.scatter(self.sochimi_tr,self.Hr_tot,label='Capacidad de Camas')
                plt.scatter(self.sochimi_tr,self.Vr_tot,label='Capacidad de Ventiladores')

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        for i in range(self.numescenarios):    
            plt.plot(self.t[i][:endD[i]],np.array(self.CH[i][:endD[i]])+np.array(self.H_sum[i][:endD[i]]),label='UCI/UTI Mov = '+str(self.inputarray[i][1]),color = 'red' ,linestyle = linestyle[i])
            plt.plot(self.t[i][:endD[i]],np.array(self.CV[i][:endD[i]])+np.array(self.V[i][:endD[i]]),label='VMI Mov = '+str(self.inputarray[i][1]),color = 'blue' ,linestyle = linestyle[i])
        self.plot(title='Necesidad total de Camas',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))




    # -------------------------------------------------------- #
    #                       Infectados                         #
    # -------------------------------------------------------- #

    # ------------------------------ #
    #       Infectados Activos       #
    # ------------------------------ #
    def plotinfectadosactivos(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True,minciencia = True):
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
            for i in range(self.numescenarios):
                plt.plot([], [], ' ', label='err: '+str(round(100*self.err_Iactives[i],2))+'%'+' Mov = '+str(self.inputarray[i][2]))    

        # Reales
        if self.realdata:
            if reales:
                if minciencia:
                    plt.scatter(self.tr,self.Ir,label='Infectados Activos reales')
                else:
                    plt.scatter(self.I_minciencia_r_tr,self.I_minciencia_r,label='Infectados Activos reales')


        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dotted',color = 'grey')

        # Infectados
        linestyle = ['dashed','solid','dashed','dotted','dotted']
        #linestyle = ['dashed','solid','dashed','solid','dotted']
        colors = ['red','blue','green','purple','black']
        #colors = ['lime','lime','purple','purple','black']
        for i in range(self.numescenarios):        
            plt.plot(self.t[i],self.I_act[i]/Isf,label='Infectados Mov = '+str(self.inputarray[i][2]),color = colors[i],linestyle=linestyle[i],linewidth=2)

        if days >0:
            plt.xlim(0,days)
        if ylim >0:
            plt.ylim(0,ylim)            
        self.plot(title = 'Activos',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)


    # ------------------------------------------ #
    #       Infectados Activos Desagregados      #
    # ------------------------------------------ #

    def plotinfectadosactivosdesagregados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
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
        # Error
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='err: '+str(round(100*err[i],2))+'%')    

        # Reales
        if self.realdata:
            if reales:
                plt.scatter(self.tr,self.Ir,label='Infectados Activos reales')

        # Infectados
        for i in range(self.numescenarios):        
            #plt.plot(self.t[i],I[i],label='Infectados )
            plt.plot(self.t[i],self.I_as[i],label='Activos asintomáticos Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_mi[i],label='Activos Mild Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_se[i],label='Activos Severos Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_cr[i],label='Activos Criticos Mov = '+str(self.inputarray[i][2]))        

        if days >0:
            plt.xlim(0,days)
        if ylim >0:
            plt.ylim(0,ylim)
        self.plot(title = 'Infectados Activos desagregados',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))

    # -------------------------------- #
    #       Infectados Acumulados      #
    # -------------------------------- #
    # No esta listo
    def plotinfectadosacumulados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
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
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='err: '+str(round(100*self.err[i],2))+'%')    

        # Reales
        if self.realdata:
            if reales:
                plt.scatter(self.I_ac_r_tr,self.I_ac_r,label='Infectados Activos reales')

        # Infectados
        for i in range(self.numescenarios):        
            plt.plot(self.t[i],self.Iac[i]/Isf,label='Infectados Mov = '+str(self.inputarray[i][2]))

        if days >0:
            plt.xlim(0,days)
        self.plot(title = 'Infectados Acumulados',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))  


    # -------------------------------------------- #
    #       Infectados Acumulados  Desagregados    #
    # -------------------------------------------- #
    
    def plotinfectadosacumuladosdesagregados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
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
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='err: '+str(round(100*self.err[i],2))+'%')    

        # Reales
        #if reales:
        #    plt.scatter(tr,Ir,label='Infectados Activos reales')

        # Infectados
        
        for i in range(self.numescenarios):        
            #plt.plot(self.t[i],I[i],label='Infectados )
            plt.plot(self.t[i],self.I_as_ac[i],label='Acumulados asintomáticos Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_mi_ac[i],label='Acumulados Mild Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_se_ac[i],label='Acumulados Severos Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_cr_ac[i],label='Acumulados Criticos Mov = '+str(self.inputarray[i][2]))        

        if days >0:
            plt.xlim(0,days)
        if ylim >0:
            plt.ylim(0,ylim)
        self.plot(title = 'Infectados Acumulados desagregados',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))        

    # ------------------------------ #
    #       Infectados Diarios       #
    # ------------------------------ #    
    def plotinfectadosdiarios(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
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


        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')       

        Idiarios = [self.I_as_d[i]+self.I_mi_d[i]+self.I_se_d[i]+self.I_cr_d[i] for i in range(self.numescenarios)]

       

        for i in range(self.numescenarios):        
            plt.plot(self.t[i][:(endD[i])],Idiarios[i][:endD[i]],label='Infectados Mov = '+str(self.inputarray[i][2]))

        if days >0:
            plt.xlim(0,days)
        self.plot(title = 'Infectados Diarios',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))          
        

    # ------------------------------------------ #
    #       Infectados Diarios  Desagregados     #
    # ------------------------------------------ #
    
    def plotinfectadosdiariosdesagregados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
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
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='err: '+str(round(100*self.err[i],2))+'%')    

        # Reales
        #if reales:
        #    plt.scatter(tr,Ir,label='Infectados Activos reales')

        # Infectados
        
        for i in range(self.numescenarios):        
            #plt.plot(self.t[i],I[i],label='Infectados )
            plt.plot(self.t[i],self.I_as_d[i],label='Diarios asintomáticos Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_mi_d[i],label='Diarios Mild Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_se_d[i],label='Diarios Severos Mov = '+str(self.inputarray[i][2]))
            plt.plot(self.t[i],self.I_cr_d[i],label='Diarios Criticos Mov = '+str(self.inputarray[i][2]))        

        if days >0:
            plt.xlim(0,days)
        if ylim >0:
            plt.ylim(0,ylim)
        self.plot(title = 'Infectados Diarios desagregados',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))        



    # ------------------------------------------------------------------------------------------------------- #
    #                                            Fallecidos                                                   #
    # ------------------------------------------------------------------------------------------------------- #

    # --------------------------------- #
    #      Fallecidos  acumulados       #
    # --------------------------------- #
    def plotfallecidosacumulados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True):
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

      
        # Inicio cuarentena general
        plt.axvline(x=self.inputarray[0][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dashed',color = 'grey')
        
        # ----------- #
        #     Plot    #
        # ----------- #
        # Parametros 
        plt.plot([], [], ' ', label='beta: '+str(self.beta))
        plt.plot([], [], ' ', label='mu: '+str(self.mu))
        plt.plot([], [], ' ', label='factor de escala: '+str(self.ScaleFactor))

        # Fecha de Peak
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))        
        # Error
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'err: '+str(round(100*err[i],2))+'%')

        # Datos reales
        if self.realdata:
            if reales:
                plt.scatter(self.Br_tr,self.Br,label='Fallecidos reales')
                #plt.scatter(self.ED_tr,self.ED_RM_ac,label='Fallecidos excesivos proyectados')

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        #linestyle = ['dashed','solid','dashed','solid','dotted']
        colors = ['red','blue','green','purple','black']
        #colors = ['lime','lime','purple','purple','black']
        for i in range(self.numescenarios):
            plt.plot(self.t[i][:endD[i]],self.B[i][:endD[i]]/norm,label='Fallecidos Mov = '+str(self.inputarray[i][2]),color = colors[i],linestyle=linestyle[i],linewidth=2)

        plt.xlim(0,days)   
        if ylim >0:
            plt.ylim(0,ylim)

        self.plot(title = 'Fallecidos',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)
        


    # ----------------------------- #
    #       Fallecidos diarios      #
    # ----------------------------- #
    def plotfallecidosdiarios(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
    

        Isf = 1    
        if scalefactor:
            Isf = self.ScaleFactor

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        for i in range(self.numescenarios):
            plt.plot(self.t[i],self.D[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])
        self.plot(title = 'Fallecidos diarios',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))

    # ------------------------------------------- #
    #       Fallecidos Desagregados Acumulados    #
    # ------------------------------------------- #
    def plotfallecidosdesagregados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True, accumulated = True):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
    

        Isf = 1    
        if scalefactor:
            Isf = self.ScaleFactor

        linestyle = ['solid','dashed','solid','dashed','dotted','dotted']
        if accumulated:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.H_crD[i]/Isf,label='H_crD - Mov = '+str(self.inputarray[i][2]),linestyle = linestyle[i])
                plt.plot(self.t[i],self.I_crD[i]/Isf,label='I_crD - Mov = '+str(self.inputarray[i][2]),linestyle = linestyle[i])
                plt.plot(self.t[i],self.I_seD[i]/Isf,label='I_seD - Mov = '+str(self.inputarray[i][2]),linestyle = linestyle[i])
                plt.plot(self.t[i],self.VD[i]/Isf,label='VD - Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])
        else:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.H_crD_d[i]/Isf,label='H_crD - Mov = '+str(self.inputarray[i][2]),linestyle = linestyle[i])
                plt.plot(self.t[i],self.I_crD_d[i]/Isf,label='I_crD - Mov = '+str(self.inputarray[i][2]),linestyle = linestyle[i])
                plt.plot(self.t[i],self.I_seD_d[i]/Isf,label='I_seD - Mov = '+str(self.inputarray[i][2]),linestyle = linestyle[i])
                plt.plot(self.t[i],self.VD_d[i]/Isf,label='VD - Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])        
        self.plot(title = 'Hospitalized Deaths',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))


    # ------------------------------------------------ #
    #     Infectados Críticos Fallecidos acumulados    #
    # ------------------------------------------------ #
    def plotfallecidosIcriticos(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True, accumulated = True):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
    

        Isf = 1    
        if scalefactor:
            Isf = self.ScaleFactor

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        if accumulated:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.I_crD[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])
        else:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.I_crD_d[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])                
        self.plot(title = 'Critical Infected Deaths',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))

    # ------------------------------------------------ #
    #     Infectados Severos Fallecidos acumulados    #
    # ------------------------------------------------ #
    def plotfallecidosIseveros(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True, accumulated = True):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
    

        Isf = 1    
        if scalefactor:
            Isf = self.ScaleFactor

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        if accumulated:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.I_seD[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])
        else:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.I_seD_d[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])                            
        self.plot(title = 'Severe Infected Deaths',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))


    # --------------------------------------- #
    #     Ventilados Fallecidos acumulados    #
    # --------------------------------------- #
    def plotfallecidosventilados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True, accumulated = True):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
    

        Isf = 1    
        if scalefactor:
            Isf = self.ScaleFactor

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        if accumulated:        
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.VD[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])
        else:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.VD_d[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])            
        self.plot(title = 'Ventilated Deaths',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))

    # ------------------------------------------- #
    #     Hospitalizados Fallecidos acumulados    #
    # ------------------------------------------- #
    def plotfallecidoshospitalizados(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True, accumulated = True):
        # -------- #
        #   Time   #
        # -------- #
        if days == 0:
            days = (enddate-self.initdate).days      
        elif days < 0:
            days = self.tsim     
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
    

        Isf = 1    
        if scalefactor:
            Isf = self.ScaleFactor

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        if accumulated:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.H_crD[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])
        else:
            for i in range(self.numescenarios):
                plt.plot(self.t[i],self.H_crD_d[i]/Isf,label='Mov = '+str(self.inputarray[i][2]),color = 'black' ,linestyle = linestyle[i])            
        self.plot(title = 'Hospitalized Deaths',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))


    # ---------------------- #
    #       Letalidad        #
    # ---------------------- #

    def plotletalidad(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True):
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
        plt.plot([], [], ' ', label='beta: '+str(self.beta))
        plt.plot([], [], ' ', label='mu: '+str(self.mu))
        plt.plot([], [], ' ', label='factor de escala: '+str(self.ScaleFactor))
        
        # Inicio cuarentena general
        for i in range(self.numescenarios):
            plt.axvline(x=self.inputarray[i][4],linestyle = 'dashed',color = 'grey')
        # Fin cuarentena general
        plt.axvline(x=self.inputarray[0][5],linestyle = 'dotted',color = 'grey')

        # Infectados
        linestyle = ['dashed','solid','dashed','dotted','dotted']
        #linestyle = ['dashed','solid','dashed','solid']
        colors = ['red','blue','green','purple','black']
        #colors = ['lime','lime','purple','purple','black']

        # Fecha de Peak
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        
        # Error
        #for i in range(self.numescenarios):
        #    plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'err: '+str(round(100*err[i],2))+'%')


        # Datos reales
        #if reales:
        #    plt.scatter(Br_tr,Br,label='Fallecidos reales')
        #    plt.scatter(self.ED_tr,self.ED_RM_ac,label='Fallecidos excesivos proyectados')

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        linestyle = ['dashed','solid','dashed','solid','dotted']
        colors = ['red','blue','green','purple','black']
        colors = ['lime','lime','purple','purple','black']
        for i in range(self.numescenarios):
            plt.plot(self.t[i],100*self.B[i]/self.Iac[i],label='Mov=['+str(self.inputarray[i][2])+','+str(self.inputarray[i][1])+']' ,color=colors[i],linestyle=linestyle[i],linewidth=2)
            
        plt.xlim(0,days)   
        if ylim >0:
            plt.ylim(0,ylim)

        self.plot(title = 'Letalidad',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))
        




    # ------------------ #
    #     Expuestos      #
    # ------------------ #
    def plotexpuestos(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False,legend=True):
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
        plt.plot([], [], ' ', label='beta: '+str(self.beta))
        plt.plot([], [], ' ', label='mu: '+str(self.mu))
        plt.plot([], [], ' ', label='factor de escala: '+str(self.ScaleFactor))

        # Fecha de Peak
        for i in range(self.numescenarios):
            plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        

        linestyle = ['dashed','solid','dashed','dotted','dotted']
        for i in range(self.numescenarios):
            plt.plot(self.t[i][:endD[i]],self.E[i][:endD[i]],label='Expuestos Mov = '+str(self.inputarray[i][2]),color = 'blue',linestyle=linestyle[i])
            plt.plot(self.t[i][:endD[i]],self.E_sy[i][:endD[i]],label='Expuestos sintomáticos Mov = '+str(self.inputarray[i][2]),color = 'red',linestyle=linestyle[i])
            
        plt.xlim(0,days)   
        if ylim >0:
            plt.ylim(0,ylim)

        self.plot(title = 'Expuestos',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'),legend=legend)
        


    # -------------------- #
    #     Curvas SEIR      #
    # -------------------- #
    def plotseird(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,seird = [1,1,1,1,1]):
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
        plt.plot([], [], ' ', label='beta: '+str(self.beta))
        plt.plot([], [], ' ', label='mu: '+str(self.mu))
        plt.plot([], [], ' ', label='factor de escala: '+str(self.ScaleFactor))

        # Fecha de Peak
        for i in range(self.numescenarios):
            plt.plot([], [], ' ', label='Mov='+str(self.inputarray[i][2])+'Peak='+self.peak_date[i].strftime('%Y-%m-%d'))
        

        linestyle = ['solid','dashed','dotted','solid','dashed','dotted']
        for i in range(self.numescenarios):        
            plt.plot(self.t[i],self.S[i],label='Susceptibles Mov = '+str(self.inputarray[i][2]),linestyle=linestyle[i])
            plt.plot(self.t[i],self.I[i],label='Infectados Mov = '+str(self.inputarray[i][2]),linestyle=linestyle[i])        
            plt.plot(self.t[i],self.E[i],label='Expuestos Mov = '+str(self.inputarray[i][2]),linestyle=linestyle[i])
            plt.plot(self.t[i],self.R[i],label='Recuperados Mov = '+str(self.inputarray[i][2]),linestyle=linestyle[i])
            #plt.plot(self.t[i],D[i],label='Muertos diarios Mov = '+str(inputarray[i][2]),linestyle=linestyle[i])
            plt.plot(self.t[i],self.B[i],label='Enterrados Mov = '+str(self.inputarray[i][2]),linestyle=linestyle[i])
            
        plt.xlim(0,days)   
        if ylim >0:
            plt.ylim(0,ylim)

        self.plot(title = 'Curvas SEIR',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))
        

    # ------------------------- #
    #       Plot Cuarentenas    #
    # ------------------------- #
    def plotcuarentenas(self,enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False):
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
        for i in range(self.numescenarios):        
            plt.plot(self.t[i],self.quarantines[i],label='Cuarentena Mov = '+str(self.inputarray[i][2]))

        if days >0:
            plt.xlim(0,days)
        self.plot(title = 'Cuarentenas',xlabel='Dias desde '+datetime.strftime(self.initdate,'%Y-%m-%d'))  


    """
    # ------------------------------------------ #
    #       Graficos para parametrización        #
    # ------------------------------------------ #
    
    # ----------------------------------------- #
    #       Curvas Expuestos/Infectados         #
    # ----------------------------------------- #
    def plotexpuestosinfectados(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,seird = [1,1,1,1,1]):
        enddate =  datetime(2020,6,30)
        days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]  
        EIrate = [self.E[i]/self.I_sum[i] for i in range(self.numescenarios)]

        for i in range(self.numescenarios):    
            plt.plot(self.t[i][:endD[i]],EIrate[:endD[i]],label='Tasa Expuestos/Infectados')
        self.plot(title='Expuestos/infectados - mu ='+str(self.mu)+' beta='+str(self.beta))




    # ------------------------ #
    #       Curvas H/I         #
    # ------------------------ #
    def plothospitalizadosinfectados(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,seird = [1,1,1,1,1]):
        #initday = self.initdate#date(2020,3,15)
        enddate =  datetime(2020,6,30)
        days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]

        HIrate = [self.H_sum[i]/self.I_sum[i] for i in range(self.numescenarios)]
        for i in range(self.numescenarios):    
            plt.plot(self.t[i][:endD[i]],HIrate[:endD[i]],label='Tasa Expuestos/Infectados')
        self.plot(title='H/I - mu ='+str(self.mu)+' beta='+str(self.beta))


    # ------------------------ #
    #       Curvas V/I         #
    # ------------------------ #
    def plotventiladosinfectados(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,seird = [1,1,1,1,1]):
        #initday = self.initdate#date(2020,3,15)
        enddate =  datetime(2020,6,30)
        days = (enddate-self.initdate).days
        days = 200
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        
        VIrate = [self.V[i]/self.I_sum[i] for i in range(self.numescenarios)]
        for i in range(self.numescenarios):     
            plt.plot(self.t[i][:endD[i]],VIrate[:endD[i]],label='Tasa Expuestos/Infectados')
        self.plot(title='V/I - mu ='+str(self.mu)+' beta='+str(self.beta))


    # ------------------------ #
    #       Curvas V/H         #
    # ------------------------- #
    def plotventiladoshospitalizados(self,enddate =  datetime(2020,7,30),days=-1, reales= True,ylim = 0,norm=1,scalefactor = False,seird = [1,1,1,1,1]):
        initday = self.initdate#date(2020,3,15)
        enddate =  datetime(2020,6,30)
        days = (enddate-self.initdate).days
        endD = [np.searchsorted(self.t[i],days) for i in range(self.numescenarios)]
        VHrate = [self.V[i]/self.H_sum[i] for i in range(self.numescenarios)]
        for i in range(self.numescenarios):        
            plt.plot(self.t[i][:endD[i]],VHrate[:endD[i]],label='Tasa Expuestos/Infectados')
        self.plot(title='V/H - mu ='+str(self.mu)+' beta='+str(self.beta))


    """