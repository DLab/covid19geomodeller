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
    """
        DLab Import Data Library
        input: 
            - tstate: [string or string list] CUT por comuna o región
            - initdate: [datetime object] Initial date
        output:
            - ImportData Object.

        Usage Example: 
        a = ImportData(tstate='13', initdate=datetime(2020,5,15))
        a.importdata()
    """    
    def __init__(self,tstate,initdate):
        self.tstate = tstate
        self.initdate = initdate
    # ------------------------------- #
    #        Importar Data Real       #
    # ------------------------------- #

    # -------------------------------- #
    #            Poblacion             #
    # -------------------------------- #

    def importPopulation(self=None,endpoint = '',tstate = ''):     
        """
            Import Population
            This Function imports the selected area population. 

            input: 
                - tstate: [string or string list] CUT por comuna o región
                - endpoint [string](opcional)
            output:
                - population [int] 

            object variables:
                self.population
                self.population_male
                self.population_female

            Example usage as function:
            population = importPopulation(tstate = '13101',initdate=datetime(2020,5,15))                
            
        """
        print('Importing Population') 

        if self:
            tstate = self.tstate
        else:
            if not tstate:
                raise Exception("State code missing")

        endpointComunas = 'http://192.168.2.223:5006/getComunas'
        endpointRegiones = 'http://192.168.2.223:5006/getStates'
        endpointSS = 'http://192.168.2.223:5006/getHealthServices'
        
        county = pd.read_json(endpointComunas)
        #regiones = pd.read_json(endpointRegiones)
        #ServicioSalud =  pd.read_json(endpointSS)

        if type(tstate) == list:
            population = 0
            population_male = 0
            population_female = 0
            for i in tstate:
                if len(i)==2:                    
                    population += int(county.loc[county['state'] == int(i)]['total_pop'].sum())
                    population_male += int(county.loc[county['state'] == int(i)]['male_pop'].sum())
                    population_female += int(county.loc[county['state'] == int(i)]['female_pop'].sum())
                if len(i)>2:
                    population += int(county.loc[county['county'] == int(i)]['total_pop'].sum())
                    population_male += int(county.loc[county['county'] == int(i)]['male_pop'].sum())
                    population_female += int(county.loc[county['county'] == int(i)]['female_pop'].sum())         
        else:
            if len(tstate)==2:
                population = int(county.loc[county['state'] == int(tstate)]['total_pop'].sum())
                population_male = int(county.loc[county['state'] == int(tstate)]['male_pop'].sum())
                population_female = int(county.loc[county['state'] == int(tstate)]['female_pop'].sum())
            if len(tstate)>2:
                population = int(county.loc[county['county'] == int(tstate)]['total_pop'].sum())
                population_male += int(county.loc[county['county'] == int(tstate)]['male_pop'].sum())
                population_female += int(county.loc[county['county'] == int(tstate)]['female_pop'].sum())                        
               
        if self:
            self.population = population
            self.population_male = population_male
            self.population_female = population_female
            return
        else:        
            return population

    def importPopulationMinCiencia(self=None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv',tstate = ''):     
        """
            Import Population
            This Function imports the selected area population from the Minciencia Github Repository 

            input: 
                - tstate: [string or string list] CUT por comuna o región
                - endpoint [string](opcional)
            output:
                - population [int] 

            Variables when used as Object:
                self.population

            Usage as function:
            population = importPopulationMinCiencia(tstate = '13101',initdate=datetime(2020,5,15))

        """            
        print('Importing Population') 

        if self:
            tstate = self.tstate
        else:
            if not tstate:
                raise Exception("State code missing")

        aux = pd.read_csv(endpoint)            
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
               
        if self:
            self.population = population
            return
        else:        
            return population

    # -------------------------------- #
    #          Active Infected         #
    # -------------------------------- #

    def importActiveInfected(self=None,tstate = '',initdate=None,endpoint = 'http://192.168.2.223:5006/getActiveCasesAllComunas'):
        """
            Import Active infected
            This function imports the active infected calculated as the infected that are within their 14 days since exposure.
            
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - endpoint (opcional):
                
            output: 
                - Ir: list of Active infected 
                - tr: days from simulation beginning
                - Ir_dates: dates [datetime object]

            Variables when used as Object:
                self.Ir
                self.tr
                self.Ir_dates

            Usage as function example:
                Ir,tr,Ir_dates = importActiveInfected(tstate = '13101',initdate=datetime(2020,5,15))
                
        """ 
        print('Importing Active infected')       
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
        data = pd.DataFrame(requests.get(endpoint).json()['data'])

        if type(tstate) == list:
            counties = [i for i in tstate if len(i)>2 ]
            states = [i for i in tstate if len(i)==2 ]            
            aux = []
            for i in states:
                aux.append(data.filter(regex='^'+i,axis=1))
            
            aux.append(data[counties])
            Ir = np.array(pd.concat(aux, axis=1).sum(axis=1))
        
        else:
            if len(tstate) == 2:
                Ir = np.array(data.filter(regex='^'+tstate,axis=1).sum(axis=1))
            elif len(tstate) > 2:
                Ir = np.array(data[tstate])


        dates = list(requests.get(endpoint).json()['dates'])
        Ir_dates = [datetime.strptime(dates[i][:10],'%Y-%m-%d') for i in range(len(dates))]
        index = np.where(np.array(Ir_dates) >= initdate)[0][0]     
        Ir=Ir[index:]
        Ir_dates=Ir_dates[index:]
        tr = [(Ir_dates[i]-initdate).days for i in range(len(Ir))]
        
        if self:
            self.Ir = Ir
            self.Ir_dates = Ir_dates
            self.tr = tr            
            return
        else:        
            return Ir,tr,Ir_dates

    #def importActiveInfectedMinciencia(self=None,tstate = '',initdate=None,endpoint = ""):
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
        print('Importing Active Infected by Minciencia')
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
        
        
        if self:
            self.I_minciencia_r = I_minciencia_r
            self.I_minciencia_r_dates = I_minciencia_r_dates
            self.I_minciencia_r_tr = I_minciencia_r_tr
            return
        else:        
            return I_minciencia_r, I_minciencia_r_tr, I_minciencia_r_dates 


    # -------------------------------- #
    #      Accumulated Infected        #
    # -------------------------------- #
    def importAccumulatedInfected(self=None,tstate = '',initdate = None, endpoint = ''):     
        """
            Import acumulated infected
            This Function imports the selected area accumulated infected from CyV Endpoint
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (opcional): 
            output: 
                - I_ac_r: Real Acumulated infected
                - I_ac_r_tr: days from simulation first day
                - I_ac_r_dates: data dates

            Variables when used as Object:
                self.I_ac_r
                self.I_ac_r_tr
                self.I_ac_r_dates

            Usage as function:
                I_ac_r, I_ac_r_tr,I_ac_r_dates= importAccumulatedInfected(tstate = '13101',initdate=datetime(2020,5,15))

                
        """
        print('Importing Accumulated Infected') 
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")
        
        
        endPointTotal = 'http://192.168.2.223:5006/getTotalCasesAllComunas'
        data = pd.DataFrame(requests.get(endPointTotal).json()['data'])


        if type(tstate) == list:
            counties = [i for i in tstate if len(i)>2 ]
            states = [i for i in tstate if len(i)==2 ]            
            aux = []
            for i in states:
                aux.append(data.filter(regex='^'+i,axis=1))
            
            aux.append(data[counties])
            I_ac_r = np.array(pd.concat(aux, axis=1).sum(axis=1))
        
        else:
            if len(tstate) == 2:
                I_ac_r = np.array(data.filter(regex='^'+tstate,axis=1).sum(axis=1))
            elif len(tstate) > 2:
                I_ac_r = np.array(data[tstate])

        # Get and filter by dates
        dates = pd.DataFrame(requests.get(endPointTotal).json()['dates'])[0]
        I_ac_r_dates = [datetime.strptime(i[:10],'%Y-%m-%d') for i in dates]        
        index = np.where(np.array(I_ac_r_dates) >= initdate)[0][0] 
        I_ac_r = I_ac_r[index:]
        I_ac_r_dates = I_ac_r_dates[index:]
        I_ac_r_tr = [(I_ac_r_dates[i]-initdate).days for i in range(len(I_ac_r))]        
               
        if self:
            self.I_ac_r = I_ac_r
            self.I_ac_r_dates = I_ac_r_dates
            self.I_ac_r_tr = I_ac_r_tr
            return
        else:        
            return I_ac_r,I_ac_r_tr,I_ac_r_dates


    def importAccumulatedInfectedMinCiencia(self=None,tstate = '',initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv'):     
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
        print('Importing Accumulated Infected') 
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
               
        if self:
            self.I_ac_r = I_ac_r
            self.I_ac_r_dates = I_ac_r_dates
            self.I_ac_r_tr = I_ac_r_tr
            return
        else:        
            return I_ac_r,I_ac_r_tr,I_ac_r_dates

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
            usage 
                I_d_r, I_d_r_tr, I_d_r_dates = importDailyInfected(tstate = '13101',initdate = datetime(2020,5,15))      
        """
        print('Importing Daily Infected')

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

        
        if self:
            self.I_d_r = np.array(I_d_r_smooth)
            self.I_d_r_tr = I_d_r_tr
            self.I_d_r_dates = I_d_r_dates
            self.I_d_r_raw = I_d_r
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
        print('Importing Daily infected with backpropagated correction')
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
        Import SOCHIMI data per state.
        Currently it just supports states, but soon I'll add Health Services as the minimum territorial data.
        input:
            - tstate: región
            - initdate: Fecha de inicio
            - endpoint: Data endpoint
        output:
            - sochimi, sochimi_dates, sochimi_tr, Hr, Hr_tot, Vr, Vr_tot, 
             (data frame, fechas, dias desde inicio sim, Hospitalizados, capacidad hospitalizados, ventilados, capacidad ventilados) 
        
        Normal Usage:
            sochimi, sochimi_dates, sochimi_tr, Hr, Hr_tot, Vr, Vr_tot = importSOCHIMI(tstate = '13', initdate = datetime(2020,5,15))
        """
        print('Importing Sochimi Data') 
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
        Vr_conf =  sochimi['Vmi covid19 confirmados']
        Vr_susp =  sochimi['Vmi covid19 sospechosos']
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
               

        if self:
            self.sochimi = sochimi
            self.Hr = Hr
            self.Vr = Vr
            self.Vr_conf = Vr_conf
            self.Vr_susp = Vr_susp
            self.Vr_tot = Vr_tot
            self.Hr_tot = Hr_tot
            self.sochimi_dates = sochimi_dates
            self.sochimi_tr = sochimi_tr
            return
        else:        
            return sochimi,sochimi_dates ,sochimi_tr,Hr, Hr_tot, Vr ,Vr_tot



    # ------------------ #
    #    Datos Sochimi   #
    # ------------------ #
    def importSOCHIMIMinCiencia(self=None,tstate = '', initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto48/SOCHIMI.csv'):
        """
        Import SOCHIMI data per state.
        Currently it just supports states, but soon I'll add Health Services as the minimum territorial data.
        input:
            - tstate: región
            - initdate: Fecha de inicio
            - endpoint: Data endpoint
        output:
            - sochimi, sochimi_dates, sochimi_tr, Hr, Hr_tot, Vr, Vr_tot, 
             (data frame, fechas, dias desde inicio sim, Hospitalizados, capacidad hospitalizados, ventilados, capacidad ventilados) 
        
        Normal Usage:
            sochimi, sochimi_dates, sochimi_tr, Hr, Hr_tot, Vr, Vr_tot = importSOCHIMI(tstate = '13', initdate = datetime(2020,5,15))
        """
        print('Importing Sochimi Data 2') 
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
        
        sochimi = pd.read_csv(endpoint)

        sochimi = sochimi.loc[sochimi['Codigo region'] == int(tstate)]

        UTI = sochimi.loc[sochimi['Serie'] == 'Camas ocupadas intermedio'].iloc[:,4:].sum()
        UTI_tot = sochimi.loc[sochimi['Serie'] == 'Camas totales intermedio'].iloc[:,4:].sum()
        UCI = sochimi.loc[sochimi['Serie'] == 'Camas ocupadas intensivo'].iloc[:,4:].sum()
        UCI_tot = sochimi.loc[sochimi['Serie'] == 'Camas totales intensivo'].iloc[:,4:].sum()
        
        VMI = sochimi.loc[sochimi['Serie'] == 'Vmi ocupados'].iloc[:,4:].sum()
        VMI_tot = sochimi.loc[sochimi['Serie'] == 'Vmi totales'].iloc[:,4:].sum()

        VMI_sospechoso = sochimi.loc[sochimi['Serie'] == 'Vmi covid19 sospechosos'].iloc[:,4:].sum()
        VMI_confirmado = sochimi.loc[sochimi['Serie'] == 'Vmi covid19 confirmados'].iloc[:,4:].sum()

        
        sochimi_dates = [datetime.strptime( sochimi.columns[4:].tolist()[i],'%Y-%m-%d') for i in range(len(sochimi.columns[4:].tolist()))]        
        index = np.where(np.array(sochimi_dates) >= initdate)[0][0] 

        
        UCI = list(UCI[index:])
        UCI_tot = list(UCI_tot[index:])
        UTI = list(UTI[index:])
        UTI_tot = list(UTI_tot[index:])

        VMI = list(VMI[index:])
        VMI_tot = list(VMI_tot[index:])
        VMI_sospechoso = list(VMI_sospechoso[index:])
        VMI_confirmado = list(VMI_confirmado[index:])

        sochimi_dates = list(sochimi_dates[index:])
        sochimi_tr = [(sochimi_dates[i]-initdate).days for i in range(len(sochimi_dates))]
        sochimi = sochimi[index:] 

        Hr = np.array(UCI) #+ np.array(UTI)
        Hr_tot = np.array(UCI_tot)# + np.array(UTI_tot)
        
               
        if self:
            self.UCI = UCI
            self.UCI_tot = UCI_tot
            self.UTI = UTI
            self.UTI_tot = UTI_tot
            self.VMI = VMI
            self.VMI_tot = VMI_tot
            self.VMI_sospechoso = VMI_sospechoso
            self.VMI_confirmado = VMI_confirmado
            self.sochimi_dates = sochimi_dates
            self.sochimi_tr = sochimi_tr
            self.sochimi = sochimi
            self.Hr = Hr
            self.Hr_tot = Hr_tot
            self.Vr = VMI
            self.Vr_tot = VMI_tot

            return
        else:        
            return UCI, UCI_tot, UTI, UTI_tot, VMI, VMI_tot, VMI_sospechoso, VMI_confirmado, sochimi_dates, sochimi_tr, sochimi, Hr, Hr_tot



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
            Usage:
                Br,Br_tr,Br_dates = importAcumulatedDeaths(self=None,tstate = '13',initdate = datetime(2020,5,15))
        """
        print('Importing Accumulated Deaths')

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
        
        if self:
            self.Br = Br
            self.Br_dates = Br_dates
            self.Br_tr = Br_tr
            return
        else:        
            return Br,Br_tr,Br_dates

    #self.importfallecidosacumulados = self.importAcumulatedDeaths

    # ------------------------------------ #
    #    Datos Fallecidos Hospitalizados   #
    # ------------------------------------ #
    def importDeathsHospitalized(self=None,tstate = '',initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto57/fallecidos_hospitalizados.csv' ):     
        """
            Import Acumulated Deaths
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (optional): 
            output: 
                - Dr_hosp: Hospitalized deaths
                - Dr_Nonhosp: Non Hospitalized deaths                
                - Br_hosp: Hospitalized Accumulated deaths
                - Br_Nonhosp: Non Hospitalized Accumulated deaths
                - Br_hosp_tr: Days since initdate
                - Dr_hosp_dates: Dates
            Usage:
                Dr_hosp, Dr_Nonhosp,Br_hosp, Br_Nonhosp, hosp_tr, hosp_dates =importAcumulatedDeathsHospitalized(tstate = '13',initdate = datetime(2020,5,15))
        """
        
        print('Importing Hospitalized/NonHospitalized Deaths')

        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")      
        
        endpoint = 'http://192.168.2.223:5006/getDeathsOriginAllStates'
        aux = pd.DataFrame( requests.get(endpoint).json()['data'])
       
        
        if type(tstate) == list:
            aux = aux[tstate[0][:2]]
            if len(tstate[0])>2:
                print('No information for counties, importing state information')
        else: 
            aux = aux[tstate[:2]]
            if len(tstate)>2:
                print('No information for counties, importing state information')

        
        Dr_hosp = aux.loc['fromHospital']
        Dr_Nonhosp = aux.loc['notFromHospital']
        hosp_dates = [datetime.strptime(aux.loc['dates'][i][:10], '%Y-%m-%d') for i in range(len(Dr_hosp))]
        Br_hosp = np.cumsum(Dr_hosp)
        Br_Nonhosp = np.cumsum(Dr_Nonhosp)
        
        
        index = np.where(np.array(hosp_dates) >= initdate)[0][0]
        
        Dr_hosp = Dr_hosp[index:]
        Dr_Nonhosp = Dr_Nonhosp[index:]
        
        Br_hosp = Br_hosp[index:]
        Br_Nonhosp = Br_Nonhosp[index:]

        hosp_dates = hosp_dates[index:]
        hosp_tr = [(hosp_dates[i]-initdate).days for i in range(len(hosp_dates))]


        #endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto57/fallecidos_hospitalizados.csv'
        #aux = pd.read_csv(endpoint)
        #aux = aux.loc[aux['Region'] == 'Metropolitana']
        #Dr_hosp = aux.loc[aux['Hospitalizacion'] == 'VERDADERO']['2020-09-07']
        #Dr_Nonhosp = aux.loc[aux['Hospitalizacion'] == 'FALSO']['2020-09-07']
        #hosp_dates =  [datetime.strptime(aux.loc[aux['Hospitalizacion'] == 'VERDADERO']['Fecha'].tolist()[i], '%Y-%m-%d') for i in range(len(Dr_hosp))]
        
        #Br_hosp = Dr_hosp.cumsum()
        #Br_Nonhosp = Dr_Nonhosp.cumsum()

        #index = np.where(np.array(hosp_dates) >= initdate)[0][0]
        #
        #Dr_hosp = Br_hosp[index:]
        #Dr_Nonhosp = Br_Nonhosp[index:]
        #
        #Br_hosp = Br_hosp[index:]
        #Br_Nonhosp = Br_Nonhosp[index:]

        #hosp_dates = hosp_dates[index:]
        #hosp_tr = [(hosp_dates[i]-initdate).days for i in range(len(hosp_dates))]
        #
        ##V_F = [Br_Nonhosp.iloc[i]/Br_hosp.iloc[i] for i in range(len(Br_Nonhosp))]
        
      
        if self:
            self.Dr_hosp = Dr_hosp
            self.Dr_Nonhosp = Dr_Nonhosp            
            self.Br_hosp = Br_hosp
            self.Br_Nonhosp = Br_Nonhosp

            self.hosp_dates = hosp_dates
            self.hosp_tr = hosp_tr
            return
        else:        
            return Dr_hosp, Dr_Nonhosp,Br_hosp, Br_Nonhosp, hosp_tr, hosp_dates            



    #self.importfallecidosacumulados = self.importAcumulatedDeaths


    # ----------------------------------------- #
    #          DEIS Deaths (Not ready)          #
    # ----------------------------------------- #
    def importDeathsDEIS(self=None,tstate = '',initdate = None):     
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
            Usage:
                Br,Br_tr,Br_dates = importAcumulatedDeaths(self=None,tstate = '13',initdate = datetime(2020,5,15))
        """
        print('Importing Deaths by DEIS')

        endpointreg = 'http://192.168.2.223:5006/getDeathsByState?state='
        endpointcom = 'http://192.168.2.223:5006/getDeathsByComuna?comuna='
 
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")      


        D_r_confirmed = []
        D_r_suspected = []

        if type(tstate) == list:
            if len(tstate[0])>2:
                aux = pd.read_json(endpointcom+tstate[0])
                D_r_dates = [datetime.strptime(aux['dates'][i][:10],'%Y-%m-%d') for i in range(len(aux))]
            else:
                aux = pd.read_json(endpointreg+tstate[0])
                D_r_dates = [datetime.strptime(aux['dates'][i][:10],'%Y-%m-%d') for i in range(len(aux))]

            for i in tstate:
                if len(i)>2:
                    aux = pd.read_json(endpointcom+i)
                else:
                    aux = pd.read_json(endpointreg+i)
                if len(D_r_confirmed)>1:
                    D_r_confirmed += np.array(aux['confirmed'])
                    D_r_suspected += np.array(aux['suspected'])
                else:
                    D_r_confirmed = np.array(aux['confirmed'])
                    D_r_suspected = np.array(aux['suspected'])
        else:        
            if len(tstate)>2:
                aux = pd.read_json(endpointcom+tstate)
            else:
                aux = pd.read_json(endpointreg+tstate)

            D_r_confirmed = np.array(aux['confirmed'])
            D_r_suspected = np.array(aux['suspected'])            
            D_r_dates = [datetime.strptime(aux['dates'][i][:10],'%Y-%m-%d') for i in range(len(aux))]
        


        index = np.where(np.array(D_r_dates) >= initdate)[0][0] 
        D_r_confirmed = D_r_confirmed[index:]
        D_r_suspected = D_r_suspected[index:]
        D_r_dates = D_r_dates[index:]
        D_r_tr = [(D_r_dates[i]-initdate).days for i in range(len(D_r_dates))]

        B_r_confirmed = D_r_confirmed.cumsum()
        B_r_suspected = D_r_suspected.cumsum()

        
        if self:
            self.Dr = D_r_confirmed
            self.Dr_suspected = D_r_suspected
            self.Br = B_r_confirmed
            self.Br_suspected = B_r_suspected
            self.Br_dates = D_r_dates
            self.Br_tr = D_r_tr            
            return
        else:        
            return D_r_confirmed,D_r_suspected,B_r_confirmed,B_r_suspected,D_r_dates,D_r_tr


    #self.importinfectadosactivosminciencia = self.importActiveInfectedMinciencia

    # ---------------------------------------- #
    #       Datos Subreporte de Infectados     #
    # ---------------------------------------- #
    """
    def importInfectedSubreport(self = None,tstate = '', initdate = None,endpoint = ''):
        
            Import calculated active infected subreport 
            This Function imports the  calculated active infected subreport by the given region. 
            When working with counties, we'll assume that the subreport is homogeneous for the whole region. 
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint [deprecated]: 
            output: 
                - I_ac_r: Real Acumulated infected
                - I_ac_r_tr: days from simulation first day
                - I_ac_r_dates: data dates

            Variables when used as Object:
                self.I_ac_r
                self.I_ac_r_tr
                self.I_ac_r_dates

            Usage as function:
                I_ac_r, I_ac_r_tr,I_ac_r_dates= importAccumulatedInfected(tstate = '13101',initdate=datetime(2020,5,15))

                

        
        print('Importing Subreported Infeted')        
        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")          

        # Import Active Infected Data
        Ir_endpoint = 'http://192.168.2.223:5006/getActiveCasesAllComunas'
        Ir = pd.DataFrame(requests.get(Ir_endpoint).json()['data'])
        

        # Import Subreport proportion data
        # This is only calculated for states and then for counties 
        endpoint = 'http://192.168.2.223:5006/getActiveCasesUnderreportByState'
        SR_data = pd.DataFrame(requests.get(endpoint).json()['data'])

        if type(tstate) == list:            

            # State list
            states = [i for i in tstate if len(i)==2 ]
            
            # Adding counties' states
            counties = [i for i in tstate if len(i)>2 ]
            
            # Get Subreport data for selected states
            for i in counties:
                if i[:2] not in states:
                    states.append(i[:2])
            
            aux = []
            for i in tstate:
                if len(i)>2:
                    aux_estimate = [Ir[i][j]*SR_data[i[:2]].loc['estimate'][j] for j in range(len(SR_data[i[:2]].loc['estimate']))]
                    aux_lower = 0
                    aux_upper = 0

                aux.append()
            for i in states:
                aux.append(data.filter(regex='^'+i,axis=1))
            
            aux.append(data[counties])
            I_ac_r = pd.concat(aux, axis=1).sum(axis=1)
        
        else:
            if len(tstate) == 2:
                I_ac_r = data.filter(regex='^'+tstate,axis=1)
            elif len(tstate) > 2:
                I_ac_r = data[tstate]

        #estimate = data[

        path = "../Data/subreporte.csv"
        subreporte = pd.read_csv(path)
        subreporte = subreporte.drop(columns = ['Unnamed: 0'])
        subreporte = subreporte.loc[subreporte['cut']==int(tstate[:2])]
        subreporte_date = [datetime.strptime(i, '%Y-%m-%d') for i in subreporte['fecha']]
        index = np.where(np.array(subreporte_date) >= initdate)[0][0]
        subreporte = subreporte.iloc[index:]
        subreporte_date = subreporte_date[index:]
        subreporte_tr = [(subreporte_date[i]-initdate).days for i in range(len(subreporte))]

        # Get and filter by dates
        dates = list(requests.get(endpoint).json()['dates'])
        SR_dates = [datetime.strptime(dates[i][:10],'%Y-%m-%d') for i in range(len(dates))]        
        index = np.where(np.array(Ir_dates) >= initdate)[0][0]
        SR_dates=SR_dates[index:]
        SR_tr = [(SR_dates[i]-initdate).days for i in range(len(SR_dates))]        

        SR=SR[index:]
        
        


        # Import Active infected and calculate Real Active Infected: 
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


        if self:
            self.subreporte = subreporte
            self.subreporte_estimate = np.array(subreporte['estimate'])
            self.subreporte_date = subreporte_date
            self.subreporte_tr = subreporte_tr
            return
        else:        
            return subreporte, np.array(subreporte['estimate']), subreporte_date, subreporte_tr

    #self.importSubreporte = self.importInfectedSubreport
    """
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
        print('Importing Excesive Deaths')
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
    def importAdjacencyRegional(self,tstate= '', endpoint = 'http://192.168.2.223:5006/getRegionalAdjacencyMatrix', N = 1):
        """
        Import adjacency data 
        N is the adjacency order, with 0 meaning the immediate neighbors, 1 the neighbor's neighbors, and so on.
        """
        print('Importing Regional Adjacency')
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
        Input:
            - tstate [string or list of strings]
            - N: Neighors Order
        """
        print('Importing National Adjacency')

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

    # ---------------------- #
    #    Import Lockdowns    #
    # ---------------------- #
    def importLockdowns(self,tstate= '', endpoint = '', N = 1):
        """
        Import Lockdowns data 
        So far from local data
        
        """
        print('Importing Lockdown data')
        if self:
            tstate = self.tstate            
        else:
            if not tstate:
                raise Exception("State code missing")
        
        endpoint = '../Data/Confinamiento_COVID19.xlsx'
        Data = pd.read_excel('Confinamiento_COVID19.xlsx',header=1,index_col=0) 
        
        Data[tstate].iloc[2:]
        # fig, ax = plt.subplots()
        # ax.plot(valdivia.astype(int))   
        # ax.fmt_xdata = mdates.DateFormatter('%Y-%m-%d') 
        # plt.show() 

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

    # --------------------------- #
    #    Importar toda la data    #
    # --------------------------- #

    def importdata(self):
        print('Importing General Data')
        try:
            self.importPopulation()
        except:
            print('Dlab Endpoint Error')
            self.importPopulationMinCiencia()
        try:
            self.importActiveInfected()
        except:
            print('Dlab Endpoint Error')
            self.importActiveInfectedMinciencia()
        try:
            self.importAccumulatedInfected()
        except:
            print('Dlab Endpoint Error')            
            self.importAccumulatedInfectedMinCiencia()

        self.importDailyInfected()
        #self.importSOCHIMI()
        self.importSOCHIMIMinCiencia()
        self.importAcumulatedDeaths()
        self.importActiveInfectedMinciencia()

        self.importDeathsHospitalized()
        #self.importInfectedSubreport()
        #self.importpcrpop()        
        #self.importfallecidosexcesivos()        
        print('Done')

    def help(self):
        help= """
                Import Data Library


              """
        print(help)
#self.importfallecidosacumulados()
#self.importinfectadosactivosminciencia()
#self.importsochimi()
#self.importPopulation()
#self.importinfectadosactivos()
#self.importinfectadosdiarios()
#self.importinfectadosacumulados()        
##self.importSubreporte()        
        
        




