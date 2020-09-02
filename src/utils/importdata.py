#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import requests
import pandas as pd
from datetime import datetime
from datetime import timedelta
import numpy as np
from os import path

"""   
   SEIRHDV Import data
"""

class ImportData():
    def __init__(self,tstate,initdate):
        self.tstate = tstate
        self.initdate = initdate
    # ------------------------------- #
    #        Importar Data Real       #
    # ------------------------------- #

    # -------------------------------- #
    #            Poblacion             #
    # -------------------------------- #
    def importPopulation(self=None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv',tstate = ''):     
        """
            Import Population
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - endpoint (opcional): 
        """
        aux = pd.read_csv(endpoint)
        if self:
            tstate = self.tstate
        else:
            if not tstate:
                raise Exception("State code missing")
        if type(tstate) == list:
            population = 0
            for i in tstate:
                if len(i)==2:
                    population += int(aux.loc[aux['Codigo region']==int(i)].iloc[:,4].sum())
                if len(i)>2:
                    population += int(aux.loc[aux['Codigo comuna']==int(i)].iloc[:,4].sum())            
        else:
            if len(tstate)==2:
                population = aux.loc[aux['Codigo region']==int(tstate)].iloc[:,4].sum()
            if len(tstate)>2:
                population = int(aux.loc[aux['Codigo comuna'] == int(tstate)].iloc[:,4])
        print('Import Population')        
        if self:
            self.population = population
            return
        else:        
            return population

    def importActiveInfected(self=None,tstate = '',initdate=None,endpoint = "http://192.168.2.223:5006/getActiveNewCasesByComuna?comuna="):
        """
            Import Active infected
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - endpoint (opcional): 
        """        
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")
        # ---------------------- # 
        #   Infectados Activos   #
        # ---------------------- #

        cutlistendpoint = 'http://192.168.2.220:8080/covid19/selectComunas' 
        cutlist  = pd.read_json(cutlistendpoint)[['cut','idState']]        

        actives = []
        mydict = None
        if type(tstate) == list:
            for i in tstate:
                if len(i)<=2:
                    aux = cutlist[cutlist['idState']==int(i)]
                    for j in aux['cut']:    
                        auxendpoint = endpoint+str(j).zfill(5)
                        r = requests.get(auxendpoint) 
                        mydict = r.json()
                        actives.append(mydict['actives'])
                        #data=pd.DataFrame(mydict)
                    #Ir = (np.array(actives)).sum(axis=0)
                elif len(i)>2:
                    auxendpoint = endpoint+i
                    r = requests.get(auxendpoint) 
                    mydict = r.json()
                    actives.append(mydict['actives'])
                    #Ir = np.array(mydict['actives'])
                Ir = (np.array(actives)).sum(axis=0)
        else:
            if len(tstate)<=2:
                aux = cutlist[cutlist['idState']==int(tstate)]
                for j in aux['cut']:
                    auxendpoint = endpoint+str(j).zfill(5)
                    r = requests.get(auxendpoint) 
                    mydict = r.json()
                    actives.append(mydict['actives'])
                    
                Ir = (np.array(actives)).sum(axis=0)
            elif len(tstate)>2:
                auxendpoint = endpoint+tstate
                r = requests.get(auxendpoint) 
                mydict = r.json()
                Ir = np.array(mydict['actives'])

        Ir_dates = [datetime.strptime(mydict['dates'][i][:10],'%Y-%m-%d') for i in range(len(mydict['dates']))]
        index = np.where(np.array(Ir_dates) >= initdate)[0][0]     
        Ir=Ir[index:]
        Ir_dates=Ir_dates[index:]
        tr = [(Ir_dates[i]-initdate).days for i in range(len(Ir))]
        print('Infectados Activos')
        if self:
            self.Ir = Ir
            self.Ir_dates = Ir_dates
            self.tr = tr            
            return
        else:        
            return Ir,tr,Ir_dates

    #self.importinfectadosactivos = self.importActiveInfected

    # -------------------------------- #
    #    Datos Infectados acumulados   #
    # -------------------------------- #
    def importAcumulatedInfected(self=None,tstate = '',initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv'):     
        """
            Import acumulated infected
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (opcional): 
            output: 
                - I_ac_r: Real Acumulated infected
                - I_ac_r_tr: days from simulation first day
                - I_ac_r_dates: data dates
        """
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")
        
        aux = pd.read_csv(endpoint)
        
        if type(tstate) == list:            
            I_ac_r = aux.loc[aux['Codigo region'].isin(tstate)].iloc[:,5:-1]            
            I_ac_r = I_ac_r.append(aux.loc[aux['Codigo comuna'].isin(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        else:                        
            I_ac_r = aux.loc[aux['Codigo region']==int(tstate)].iloc[:,5:-1]
            I_ac_r = I_ac_r.append(aux.loc[aux['Codigo comuna']==int(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        
        I_ac_r_dates = [datetime.strptime(I_ac_r.index[i],'%Y-%m-%d') for i in range(len(I_ac_r))]
        index = np.where(np.array(I_ac_r_dates) >= initdate)[0][0] 
        I_ac_r = I_ac_r[index:]
        I_ac_r_dates = I_ac_r_dates[index:]
        I_ac_r_tr = [(I_ac_r_dates[i]-initdate).days for i in range(len(I_ac_r))]        
        print('Infectados Acumulados')        
        if self:
            self.I_ac_r = I_ac_r
            self.I_ac_r_dates = I_ac_r_dates
            self.I_ac_r_tr = I_ac_r_tr
            return
        else:        
            return I_ac_r,I_ac_r_tr,I_ac_r_dates

    #self.importinfectadosacumulados = self.importAcumulatedInfected


    # -------------------------------- #
    #      Daily infected Smoothed     #
    # -------------------------------- #
    # Created by Felipe Castillo
    def importDailyInfected(self=None,tstate = '',initdate = None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv' ):     
        """
            Import daily infected
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (optional): 
            output: 
                - I_d_r: Real Daily infected
                - I_d_r_tr: days from simulation first day
                - I_d_r_dates: data dates
        """
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")
        
        aux = pd.read_csv(endpoint)
        
        if type(tstate) == list:            
            I_ac_r = aux.loc[aux['Codigo region'].isin(tstate)].iloc[:,5:-1]            
            I_ac_r = I_ac_r.append(aux.loc[aux['Codigo comuna'].isin(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        else:                        
            I_ac_r = aux.loc[aux['Codigo region']==int(tstate)].iloc[:,5:-1]
            I_ac_r = I_ac_r.append(aux.loc[aux['Codigo comuna']==int(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        
        I_ac_r_dates = [datetime.strptime(I_ac_r.index[i],'%Y-%m-%d') for i in range(len(I_ac_r))]
        index = np.where(np.array(I_ac_r_dates) >= initdate)[0][0] 
        I_ac_r = I_ac_r[index:]
        I_ac_r_dates = I_ac_r_dates[index:]
        I_ac_r_tr = [(I_ac_r_dates[i]-initdate).days for i in range(len(I_ac_r))]    

        # Create daily infected from Interpolating and diferentiating the acumulated infected
        I_d_r = np.diff(np.interp(list(range(I_ac_r_tr[-1])),I_ac_r_tr,I_ac_r))
        I_d_r_tr = list(range(len(I_d_r)))
        I_d_r_dates = [initdate + timedelta(days=i) for i in range(len(I_d_r_tr))]

        outliers_init = (datetime(2020,6,15)-initdate).days
        outliers_end = (datetime(2020,6,19)-initdate).days

        I_d_r_smooth=pd.DataFrame(I_d_r)

        I_d_r_smooth[outliers_init:outliers_end] = float((I_d_r_smooth.iloc[outliers_init-2]+I_d_r_smooth.iloc[outliers_end+1])/2)
        I_d_r_smooth = I_d_r_smooth.rolling(7, win_type='gaussian', min_periods=1, center=True).mean(std=2).round()

        print('Infectados diarios')
        if self:
            self.I_d_r = I_d_r
            self.I_d_r_tr = I_d_r_tr
            self.I_d_r_dates = I_d_r_dates
            return
        else:        
            return I_d_r_smooth, I_d_r_tr, I_d_r_dates
                

    # ----------------------------------------------------- #
    #      Daily infected Smoothed with backpropagation     #
    # ----------------------------------------------------- #
    # Created by Felipe Castillo
    def importDailyInfectedBackprop(self=None,tstate = '',initdate = None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv' ):     
        """
            Import daily infected smoothed with Backpropagation
            This function smoothes the daily infected thata and distributes homogeneously the 30k excess cases added in ~2020-6-16
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (optional): 
            output: 
                - I_d_r: Real Daily infected
                - I_d_r_tr: days from simulation first day
                - I_d_r_dates: data dates
        """
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")
        
        aux = pd.read_csv(endpoint)
        
        if type(tstate) == list:            
            I_ac_r = aux.loc[aux['Codigo region'].isin(tstate)].iloc[:,5:-1]            
            I_ac_r = I_ac_r.append(aux.loc[aux['Codigo comuna'].isin(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        else:                        
            I_ac_r = aux.loc[aux['Codigo region']==int(tstate)].iloc[:,5:-1]
            I_ac_r = I_ac_r.append(aux.loc[aux['Codigo comuna']==int(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        
        I_ac_r_dates = [datetime.strptime(I_ac_r.index[i],'%Y-%m-%d') for i in range(len(I_ac_r))]
        index = np.where(np.array(I_ac_r_dates) >= initdate)[0][0] 
        I_ac_r = I_ac_r[index:]
        I_ac_r_dates = I_ac_r_dates[index:]
        I_ac_r_tr = [(I_ac_r_dates[i]-initdate).days for i in range(len(I_ac_r))]    

        # Create daily infected from Interpolating and diferentiating the acumulated infected
        I_d_r = np.diff(np.interp(list(range(I_ac_r_tr[-1])),I_ac_r_tr,I_ac_r))
        I_d_r_tr = list(range(len(I_d_r)))
        I_d_r_dates = [initdate + timedelta(days=i) for i in range(len(I_d_r_tr))]

        outliers_init = (datetime(2020,6,15)-initdate).days
        outliers_end = (datetime(2020,6,19)-initdate).days

        I_d_r_smooth=pd.DataFrame(I_d_r)

        I_d_r_smooth[outliers_init:outliers_end] = float((I_d_r_smooth.iloc[outliers_init-2]+I_d_r_smooth.iloc[outliers_end+1])/2)
        I_d_r_smooth[:outliers_end] += np.round(I_d_r_smooth[:outliers_end]*31412/np.sum(I_d_r_smooth[:outliers_end]))
        I_d_r_smooth = I_d_r_smooth.rolling(7, win_type='gaussian', min_periods=1, center=True).mean(std=2).round()

        print('Infectados diarios')
        if self:
            self.I_d_r = I_d_r
            self.I_d_r_tr = I_d_r_tr
            self.I_d_r_dates = I_d_r_dates
            return
        else:        
            return I_d_r_smooth, I_d_r_tr, I_d_r_dates
                
    

    #self.importinfectadosdiarios = self.importDailyInfected

    # ------------------ #
    #    Datos Sochimi   #
    # ------------------ #
    def importSOCHIMI(self=None,tstate = '', initdate = None, endpoint = "http://192.168.2.223:5006/getBedsAndVentilationByState?state="):
        """
        Import SOCHIMI data por región. Falta incorporar la opción de conseguir data por sector de salud.
        input:
            - tstate: región
            - initdate: Fecha de inicio
            - endpoint: Data endpoint
        output:
            - sochimi, sochimi_dates, sochimi_tr, Hr, Hr_tot, Vr, Vr_tot, 
             (data frame, fechas, dias desde inicio sim, Hospitalizados, capacidad hospitalizados, ventilados, capacidad ventilados) 
        """
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")        
        if type(tstate) == list:
            tstate = tstate[0]
        
        endpoint = endpoint+tstate[:2]
        r = requests.get(endpoint) 
        mydict = r.json()
        sochimi=pd.DataFrame(mydict)
        
        Hr = sochimi['camas_ocupadas']
        Vr =  sochimi['vmi_ocupados']
        Vr_tot =  sochimi['vmi_totales']
        Hr_tot =  sochimi['camas_totales']
        sochimi_dates = [datetime.strptime(sochimi['dates'][i][:10],'%Y-%m-%d') for i in range(len(sochimi))]

        index = np.where(np.array(sochimi_dates) >= initdate)[0][0] 
        Hr=list(Hr[index:])
        Vr=list(Vr[index:])
        Hr_tot=list(Hr_tot[index:])
        Vr_tot=(list(Vr_tot[index:]))
        sochimi_dates = sochimi_dates[index:]
        sochimi_tr = [(sochimi_dates[i]-initdate).days for i in range(len(Hr))]
        sochimi = sochimi[index:] 
        print('Sochimi')        

        if self:
            self.sochimi = sochimi
            self.Hr = Hr
            self.Vr = Vr
            self.Vr_tot = Vr_tot
            self.Hr_tot = Hr_tot
            self.sochimi_dates = sochimi_dates
            self.sochimi_tr = sochimi_tr
            return
        else:        
            return sochimi,sochimi_dates ,sochimi_tr,Hr, Hr_tot, Vr ,Vr_tot

    # -------------------------------- #
    #    Datos Fallecidos acumulados   #
    # -------------------------------- #
    def importAcumulatedDeaths(self=None,tstate = '',initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto14/FallecidosCumulativo.csv' ):     
        """
            Import Acumulated Deaths
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (optional): 
            output: 
                - Br: Real acumulated deaths
                - Br_tr: days from simulation first day
                - Br_dates: data dates
        """
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")      
                
        cut =  ['15','01','02','03','04','05','13','06','07','16','08','09','14','10','11','12','00']
        if type(tstate) == list:
            tstate = tstate[0]        
        index = cut.index(tstate[:2])
        Br = pd.read_csv(endpoint).iloc[index][1:] 
        Br_dates = [datetime.strptime(Br.index[i],'%Y-%m-%d') for i in range(len(Br))]
        index = np.where(np.array(Br_dates) >= initdate)[0][0] 
        Br = Br[index:]
        Br_dates = Br_dates[index:]
        Br_tr = [(Br_dates[i]-initdate).days for i in range(len(Br))]
        print('Fallecidos Acumulados')
        if self:
            self.Br = Br
            self.Br_dates = Br_dates
            self.Br_tr = Br_tr
            return
        else:        
            return Br,Br_tr,Br_dates

    #self.importfallecidosacumulados = self.importAcumulatedDeaths

    # ---------------------------------------- #
    #    Datos Infectados activos Minciencia   #
    # ---------------------------------------- #
    def importActiveInfectedMinciencia(self=None,tstate = '', initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto19/CasosActivosPorComuna.csv' ):     
        """
            Import Active infected Minciencia
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (optional): 
            output: 
                - I_minciencia_r: Real acumulated deaths
                - I_minciencia_r_tr: days from simulation first day
                - I_minciencia_r_dates: data dates
        """        
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")  

        aux = pd.read_csv(endpoint)

        if type(tstate) == list:            
            I_minciencia_r = aux.loc[aux['Codigo region'].isin(tstate)].iloc[:,5:-1]            
            I_minciencia_r = I_minciencia_r.append(aux.loc[aux['Codigo comuna'].isin(tstate)].iloc[:,5:-1])
            I_minciencia_r = I_minciencia_r.sum()
        else:
            I_minciencia_r = aux.loc[aux['Codigo region']==int(tstate)].iloc[:,5:-1]
            I_minciencia_r = I_minciencia_r.append(aux.loc[aux['Codigo comuna']==int(tstate)].iloc[:,5:-1])
            I_minciencia_r = I_minciencia_r.sum()


        #self.I_minciencia_r = aux.loc[aux['Codigo region']==int(self.tstate[:2])].iloc[:,5:].sum()
        

        I_minciencia_r_dates = [datetime.strptime(I_minciencia_r.index[i],'%Y-%m-%d') for i in range(len(I_minciencia_r))]
        index = np.where(np.array(I_minciencia_r_dates) >= initdate)[0][0] 
        I_minciencia_r = I_minciencia_r[index:]
        I_minciencia_r_dates = I_minciencia_r_dates[index:]
        I_minciencia_r_tr = [(I_minciencia_r_dates[i]-initdate).days for i in range(len(I_minciencia_r))]
        print('Infectados Activos Minciencia')
        
        if self:
            self.I_minciencia_r = I_minciencia_r
            self.I_minciencia_r_dates = I_minciencia_r_dates
            self.I_minciencia_r_tr = I_minciencia_r_tr
            return
        else:        
            return I_minciencia_r, I_minciencia_r_tr, I_minciencia_r_dates 

    #self.importinfectadosactivosminciencia = self.importActiveInfectedMinciencia

    # ---------------------------------------- #
    #       Datos Subreporte de Infectados     #
    # ---------------------------------------- #
    def importInfectedSubreport(self = None,tstate = '', initdate = None):
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")          
        path = "../Data/subreporte.csv"
        subreporte = pd.read_csv(path)
        subreporte = subreporte.drop(columns = ['Unnamed: 0'])
        subreporte = subreporte.loc[subreporte['cut']==int(tstate[:2])]
        subreporte_date = [datetime.strptime(i, '%Y-%m-%d') for i in subreporte['fecha']]
        index = np.where(np.array(subreporte_date) >= initdate)[0][0]
        subreporte = subreporte.iloc[index:]
        subreporte_date = subreporte_date[index:]
        subreporte_tr = [(subreporte_date[i]-initdate).days for i in range(len(subreporte))]
        print('Subreporte')
        if self:
            self.subreporte = subreporte
            self.subreporte_estimate = np.array(subreporte['estimate'])
            self.subreporte_date = subreporte_date
            self.subreporte_tr = subreporte_tr
            return
        else:        
            return subreporte, np.array(subreporte['estimate']), subreporte_date, subreporte_tr

    #self.importSubreporte = self.importInfectedSubreport

    # -------------------------- #
    #    Fallecidos excesivos    #
    # -------------------------- #
    def importfallecidosexcesivos(self,path = '/home/samuel/Documents/Dlab/data/Excess_dead_daily.csv'):
        #path = '/home/samuel/Documents/Dlab/data/Excess_dead_daily.csv'
        
        excess_dead = pd.read_csv(path)
        self.ED_RM_df = excess_dead.loc[excess_dead['Codigo region']==int(self.tstate[:2])]      
        self.ED_RM = [self.ED_RM_df['Defunciones Covid'].iloc[i] + self.ED_RM_df['Exceso de muertes media poderada'].iloc[i] for i in range(len(self.ED_RM_df))]       

        self.ED_RM_dates = [datetime.strptime(self.ED_RM_df['Fecha'].iloc[i], '%Y-%m-%d')  for i in range(len(self.ED_RM_df))]
        index = np.where(np.array(self.ED_RM_dates) >= self.initdate)[0][0]
        enddate = max(self.ED_RM_dates)
        indexend = np.where(np.array(self.ED_RM_dates) >= enddate)[0][0]
        self.ED_RM_dates = self.ED_RM_dates[index:indexend]  
        self.ED_RM = self.ED_RM[index:indexend]
        self.ED_RM_ac = np.cumsum(self.ED_RM)
        self.ED_tr = [(self.ED_RM_dates[i]-self.initdate).days for i in range(len(self.ED_RM))]
        print('Fallecidos Excesivos')
        return


    # -------------------------------- #
    #       Datos PCR y polbación      #
    # -------------------------------- #
    def importpcrpop(self,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto7/PCR.csv'):
        cut =  ['15','01','02','03','04','05','13','06','07','16','08','09','14','10','11','12','00']
        index = cut.index(self.tstate[:2])
        self.population = pd.read_csv(endpoint).iloc[index]['Poblacion'] 
        self.pcr = pd.read_csv(endpoint).iloc[index][3:] 
        self.pcr_dates = [datetime.strptime(self.pcr.index[i],'%Y-%m-%d') for i in range(len(self.Br))]
        index = np.where(np.array(self.Br_dates) >= self.initdate)[0][0] 
        self.pcr = self.pcr[index:]
        self.pcr_dates = self.pcr_dates[index:]
        self.pcr_tr = [(self.pcr_dates[i]-self.initdate).days for i in range(len(self.Br))]
        print('PCR Y poblacion')
        return


    # ----------------------------- #
    #    Import Adjacency Matrix    #
    # ----------------------------- #
    def importAdjacencyRegional(self,tstate= '', endpoint = 'http://192.168.2.223:5006/getRegionalAdjacencyMatrix', N = 0):
        """
        Import adjacency data 
        N is the adjacency order, with 0 meaning the immediate neighbors, 1 the neighbor's neighbors, and so on.
        """
        if self:
            tstate = self.tstate            
        else:
            if not tstate:
                raise Exception("State code missing")
        Data = pd.read_json(endpoint)
        if type(tstate) == list:
            adjacency = tstate
            adjacency_lvl = tstate
        else:
            adjacency = [tstate]
            adjacency_lvl = [tstate]
        adjacency_lvl.append(Data.loc[int(tstate)][0])
        adjacency.extend((Data.loc[int(tstate)][0]))
        for j in range(1,N):
            aux = [] 
            for i in adjacency_lvl[j]:
                aux.extend(Data.loc[int(i)][0])
                # Remove duplicates                
                aux = [k for k in aux if k not in adjacency]
                aux = list(set(aux))
            if len(aux)>0:                
                adjacency_lvl.append(aux)
                adjacency.extend(aux)
            else:
                break
        if self:
            self.adjacencyR = adjacency
            self.adjacencyR_lvl = adjacency_lvl
        else:        
            return adjacency, adjacency_lvl


    def importAdjacencyNational(self,tstate= '', endpoint = 'http://192.168.2.223:5006/getNationalAdjacencyMatrix', N = 1):
        """
        Import adjacency data 
        N is the adjacency order, with 0 meaning the immediate neighbors, 1 the neighbor's neighbors, and so on.
        """
        if self:
            tstate = self.tstate            
        else:
            if not tstate:
                raise Exception("State code missing")
        Data = pd.read_json(endpoint)
        if type(tstate) == list:
            adjacency = tstate
            adjacency_lvl = tstate
        else:
            adjacency = [tstate]
            adjacency_lvl = [tstate]
        adjacency_lvl.append(Data.loc[int(tstate)][0])
        adjacency.extend((Data.loc[int(tstate)][0]))
        for j in range(1,N):
            aux = [] 
            for i in adjacency_lvl[j]:
                aux.extend(Data.loc[int(i)][0])
                # Remove duplicates                
                aux = [k for k in aux if k not in adjacency]
                aux = list(set(aux))
            if len(aux)>0:                
                adjacency_lvl.append(aux)
                adjacency.extend(aux)
            else:
                break
        if self:
            self.adjacency = adjacency
            self.adjacency_lvl = adjacency_lvl
        else:        
            return adjacency, adjacency_lvl

    # --------------------------- #
    #    Importar toda la data    #
    # --------------------------- #

    def importdata(self):
        print('Importando Datos')
        self.importPopulation()
        self.importActiveInfected()
        self.importAcumulatedInfected()
        self.importDailyInfected()
        self.importSOCHIMI()
        self.importAcumulatedDeaths()
        self.importActiveInfectedMinciencia()
        #self.importInfectedSubreport()
        #self.importpcrpop()        
        #self.importfallecidosexcesivos()        
        print('Done')


#self.importfallecidosacumulados()
#self.importinfectadosactivosminciencia()
#self.importsochimi()
#self.importPopulation()
#self.importinfectadosactivos()
#self.importinfectadosdiarios()
#self.importinfectadosacumulados()        
##self.importSubreporte()        
        
        
        

    # -------------------------------- #
    #    Datos Infectados diarios      #
    # -------------------------------- #
    # Falta interpolar
    """
    def importDailyInfected(self=None,tstate = '',initdate = None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv' ):     
 
            #Import daily infected
            #input: 
            #    - tstate: [string or string list] CUT por comuna o región
            #    - initdate: datetime object with the initial date
            #    - endpoint (optional): 
            #output: 
            #    - I_d_r: Real Daily infected
            #    - I_d_r_tr: days from simulation first day
            #    - I_d_r_dates: data dates
        
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")
        
        aux = pd.read_csv(endpoint)
        
        if type(tstate) == list:            
            I_ac_r = aux.loc[aux['Codigo region'].isin(tstate)].iloc[:,5:-1]            
            I_ac_r = I_ac_r.append(aux.loc[aux['Codigo comuna'].isin(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        else:                        
            I_ac_r = aux.loc[aux['Codigo region']==int(tstate)].iloc[:,5:-1]
            I_ac_r = I_ac_r.append(aux.loc[aux['Codigo comuna']==int(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        
        I_ac_r_dates = [datetime.strptime(I_ac_r.index[i],'%Y-%m-%d') for i in range(len(I_ac_r))]
        index = np.where(np.array(I_ac_r_dates) >= initdate)[0][0] 
        I_ac_r = I_ac_r[index:]
        I_ac_r_dates = I_ac_r_dates[index:]
        I_ac_r_tr = [(I_ac_r_dates[i]-initdate).days for i in range(len(I_ac_r))]    

        I_d_r = np.diff(np.interp(list(range(I_ac_r_tr[-1])),I_ac_r_tr,I_ac_r))
        I_d_r_tr = list(range(len(I_d_r)))
        I_d_r_dates = [initdate + timedelta(days=i) for i in range(len(I_d_r_tr))]

        print('Infectados diarios')
        if self:
            self.I_d_r = I_d_r
            self.I_d_r_tr = I_d_r_tr
            self.I_d_r_dates = I_d_r_dates
            return
        else:        
            return I_d_r, I_d_r_tr, I_d_r_dates

    """
