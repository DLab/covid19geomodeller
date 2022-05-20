#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
from pandas.core.base import DataError
import requests
from requests.auth import HTTPBasicAuth
import pandas as pd
from datetime import datetime
from datetime import date
from datetime import timedelta
import numpy as np
from os import path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates 

# cv19gm libraries 
import os
path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

import sys
#from pathlib import Path
sys.path.insert(1, path)
from cv19gm.utils.cv19timeutils import timeJStoPy 

"""   
   Import and analyse data 


    To Do:
        * Save data to file 
        * Load data from file
        * Creación de df de información (self.info) sobre los datos (self.data)
        * Estandarizar nombres de variables de hospitalización
        * Expandir creación de dataframe a todas las funciones de importación
        * Función que devuelva serie de tiempo sin nulls. Puede ser necesario para funciones de análisis de datos
        * Manually add data ()

    self.resume:
        country,state,county,healthservice, population, initdate, update_date, dataindex
    self.datainfo:
        endpoint,initdate,enddate,update_date,sim_var, complete name, additional info.


    #CUT Regiones:
    
    Names = ['Arica y Parinacota','Tarapacá','Antofagasta','Atacama','Coquimbo','Valparaíso','Metropolitana','O’Higgins','Maule','Ñuble','Biobío','Araucanía','Los Ríos','Los Lagos','Aysén','Magallanes']
    cut =   ['15',                '01',      '02',         '03',     '04',      '05',        '13',           '06',       '07',   '16',   '08',    '09',       '14',      '10',       '11',   '12'] 

    nametocut = {}  
    cuttoname = {} 
    for i in range(len(Names)): 
        nametocut.update({Names[i]:cut[i]})  
        cuttoname.update({cut[i]:Names[i]})         

"""

def help():
    aux = """
            Functions: 

          """
    print(aux)


def dfappend(df,values,days,dataname):
    length = len(df)
    aux = [np.nan for i in range(length)]
    iterator = iter(values)
    for i in days:
        aux[i] = next(iterator)
    auxdf = pd.DataFrame({dataname:aux})
    return pd.concat([df,auxdf],axis=1)

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
    def __init__(self, initdate = None,enddate = datetime.now(), localdata = None):
        
        if localdata:
            if type(localdata) == str:
                self.data = pd.read_csv(localdata)
                #self.load(localdata)
            elif type(localdata) == pd.core.frame.DataFrame:
                self.data = localdata
            else:
                raise DataError
            
        else:
            tdelta = enddate - initdate 
            dates = [initdate + timedelta(days=i) for i in range(tdelta.days + 1) ] 
            days = [i for i in range(tdelta.days + 1)] 
            self.data = pd.DataFrame({'days':days,'dates':dates})

        
        self.tstate = ['15','01','02','03','04','05','13','06','07','16','08','09','14','10','11','12'] 
        self.initdate = initdate
        
    # --------------------------- #
    #    Importar toda la data    #
    # --------------------------- #

    def import_data(self):  
        self.imp_population_mcyt()        
        self.imp_infected_active_mcyt()                 
        self.imp_infected_accumulated_mcyt()
        self.imp_infected_daily_mcyt() 
        self.imp_infected_mcyt_infdiario() 
        self.imp_hosp_icu_mcyt() 
        self.imp_deaths_deis_mcyt() 
        self.imp_vacc_infected_mcyt() 
        self.imp_vaccinated_mcyt() 
        print('Done')

    def calculateIC(self,pIv_det=1,pV2=1,infected='InformeEPI',pOmicron=1):
        print('Creating Initial conditions')

        R0 = self.I_ac_r[0] - self.I_r[0] - self.Br[0] - self.UCI_use_covid[0]
        Iv0_det = self.Iv0
        Iv0 = Iv0_det/pIv_det        
        Iv_d0 = self.Iv_d[0]/pIv_det
        Iv_ac0 = self.Iv_ac[0]/pIv_det

        Iv_d0_det = self.Iv_d[0]
        Iv_ac0_det = self.Iv_ac[0]

        self.Sv0 = self.v3[0] + self.v2[0]*pV2 - Iv_ac0

        if infected == 'InformeEPI':
            I_d = self.I_d_r[0]
            I_ac = self.I_ac_r[0]

        elif infected == 'InformeDiario':
            I_d = self.I_d_infd_rmean[0]
            I_ac = self.I_ac_infd[0]            


        self.IC =  {'initdate':self.initdate.strftime("%Y-%m-%d"),'population':self.population,'Sv':self.Sv0,'I':self.I_r[0]*pOmicron,'I_d':I_d*pOmicron,
        'I_ac':I_ac,'Iv':Iv0*pOmicron,'Iv_d':Iv_d0*pOmicron,'Iv_ac':Iv_ac0,'H_cap':self.UCI_capacity[0]-self.UCI_use_noncovid[0],'H':self.UCI_use_covid[0],
        'H_d':self.UCI_d[0],'D':self.Br[0],'D_d':self.Dr[0],'R':R0}                
        print('Done')


    def calculateIC_SEIR(self,infected='InformeEPI',CrossImmunity=1):
        print('Creating Initial conditions')

        R0 = CrossImmunity*(self.I_ac_r[0] - self.I_r[0] - self.Br[0])

        if infected == 'InformeEPI':
            I_d_det = self.I_d_r[0]
            I_ac_det = self.I_ac_r[0]

        elif infected == 'InformeDiario':
            I_d_det = self.I_d_infd_rmean[0]
            I_ac_det = self.I_ac_infd[0]            


        self.IC =  {'initdate':self.initdate.strftime("%Y-%m-%d"),'population':self.population,'I_det':self.I_r[0],'I_d_det':I_d_det,'I_ac_det':I_ac_det,'R':R0}                
        print('Done')

    # ---------------------------- #
    #         Save Data            #
    # ---------------------------- #
    def save(self,filename):        
        print('Saving Data')
        self.data.to_csv(filename)
        return

    # ---------------------------- #
    #         Load Data            #
    # ---------------------------- #
    def load(self,filename):
        # Drop na values
        print('Loading Data')
        self.data = pd.read_csv(filename)
        # Process the data - create de corresponding variables and dataframes:
        return


    # ------------------------------ #
    #       Add Data Manually        #
    # ------------------------------ #
    def add_data(self,data,dates=None,days=None,initdate=None):
        print('Adding Data')
        # We use initdate in the case they give a days vector and it doesn't fit the actual data size        
        
               
        return

    # -------------------------------- #
    #            Population            #
    # -------------------------------- #

    def imp_population_mcyt(self=None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv',tstate = ''):     
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


    # ---------------------------------------- #
    #    Datos Infectados activos Minciencia   #
    # ---------------------------------------- #
    def imp_infected_active_mcyt(self=None,tstate = '', initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto19/CasosActivosPorComuna.csv',name = 'I' ):     
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
            # Update database
            print('updating database')
            self.data = dfappend(self.data,I_minciencia_r,I_minciencia_r_tr,name)            
            self.I_r = I_minciencia_r
            self.I_r_dates = I_minciencia_r_dates
            self.I_r_tr = I_minciencia_r_tr
            return
        else:        
            return I_minciencia_r, I_minciencia_r_tr, I_minciencia_r_dates 

 
    # -------------------------------- #
    #      Accumulated Infected        #
    # -------------------------------- #
    def imp_infected_accumulated_mcyt(self=None,tstate = '',initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv', name = 'I_ac'):     
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
            # Update database
            print('updating database')
            self.data = dfappend(self.data,I_ac_r,I_ac_r_tr,name)              
            self.I_ac_r = np.array(I_ac_r)
            self.I_ac_r_dates = I_ac_r_dates
            self.I_ac_r_tr = I_ac_r_tr
            return
        else:        
            return np.array(I_ac_r),I_ac_r_tr,I_ac_r_dates


    # ------------------------------------------ #
    #      Daily infected Smoothed MinCiencia    #
    # ------------------------------------------ #
    def imp_infected_daily_mcyt(self=None,tstate = '',initdate = None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv', name = 'I_d' ):     
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
                        
        data = pd.read_csv(endpoint)


        if type(tstate) == list:            
            I_ac_r = data.loc[data['Codigo region'].isin(tstate)].iloc[:,5:-1]            
            I_ac_r = I_ac_r.append(data.loc[data['Codigo comuna'].isin(tstate)].iloc[:,5:-1])
            I_ac_r = I_ac_r.sum()
        else:                        
            I_ac_r = data.loc[data['Codigo region']==int(tstate)].iloc[:,5:-1]
            I_ac_r = I_ac_r.append(data.loc[data['Codigo comuna']==int(tstate)].iloc[:,5:-1])
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

        if initdate < datetime(2020,6,15):
            outliers_init = (datetime(2020,6,15)-initdate).days
            outliers_end = (datetime(2020,6,19)-initdate).days

        I_d_r_smooth=pd.DataFrame(I_d_r)
        I_d_r_smooth = I_d_r_smooth.rolling(7, win_type='gaussian', min_periods=1, center=True).mean(std=2).round()

        
        if self:
            # Update database
            print('updating database')
            self.data = dfappend(self.data,np.array(I_d_r_smooth[0]),I_d_r_tr,name)              
            self.I_d_r = np.array(I_d_r_smooth[0])
            self.I_d_r_tr = I_d_r_tr
            self.I_d_r_dates = I_d_r_dates
            self.I_d_r_raw = I_d_r
            return
        else:        
            return np.array(I_d_r_smooth[0]), I_d_r_tr, I_d_r_dates



    # ----------------------------------------------------- #
    #            Daily infected "informe diario"            #
    # ----------------------------------------------------- #
    def imp_infected_mcyt_infdiario(self=None,initdate = None,endpoint = 'https://github.com/MinCiencia/Datos-COVID19/raw/master/output/producto5/TotalesNacionales.csv',rav=7):
        """Import national infected from the "informe diario"

        Args:
            self (object, optional): Object to add the data. Defaults to None.
            initdate (datetime.datetime, optional): initial date. Defaults to None.
            endpoint (str, optional): endpoint. Defaults to 'https://github.com/MinCiencia/Datos-COVID19/raw/master/output/producto5/TotalesNacionales.csv'.

        Raises:
            Exception: [description]
            Exception: [description]

        Returns:
            [type]: [description]
        """
        print('Importing Daily infected with backpropagated correction')
        if self:
            initdate = self.initdate
        else:
            if not initdate:
                raise Exception("Initial date missing")
        
        aux = pd.read_csv(endpoint)
        
        I_ac_infd = aux.iloc[1,1:]
        I_d_infd = aux.iloc[6,1:]
        dates_infd = [datetime.strptime(aux.columns[1:][i],'%Y-%m-%d') for i in range(len(aux.columns[1:]))]
        
        I_d_infd_rmean = I_d_infd.rolling(window=rav,win_type='gaussian', min_periods=1, center=True).mean(std=2).round()

        index = np.where(np.array(dates_infd) >= initdate)[0][0] 
        I_d_infd = I_d_infd[index:]
        I_ac_infd = I_ac_infd[index:]
        dates_infd = dates_infd[index:]
        I_d_infd_rmean = I_d_infd_rmean[index:]

        t_infd = [(dates_infd[i]-initdate).days for i in range(len(dates_infd))]


        if self:
            # Update database
            print('updating database')
            #self.data = dfappend(self.data,I_d_r_smooth,I_d_r_tr,name)              
            self.I_d_infd = I_d_infd
            self.I_ac_infd = I_ac_infd
            self.I_d_infd_rmean = I_d_infd_rmean
            self.t_infd = t_infd
            self.dates_infd = dates_infd
            return
        else:        
            return I_d_infd, I_ac_infd, I_d_infd_rmean, t_infd, dates_infd
                

    # --------------------------------------------------------- #
    #                 Ocupación Hospitalaria                    #
    # --------------------------------------------------------- #    

    # ----------------------------- #
    #    Datos Ocupacion de Camas   #
    # ----------------------------- #
    def imp_hosp_icu_mcyt(self=None,tstate = '', initdate = None, endpoint = "https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto58/Camas_UCI_diarias.csv",user=None,password=None, name = ['UCI_capacity','UCI_use_covid','UCI_use_noncovid']):
        """
        Import ICU Bed Occupation data per region.
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

        if self:
            tstate = self.tstate
            initdate = self.initdate
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")

        cuttoname = {'15': 'Arica y Parinacota', '01': 'Tarapacá', '02': 'Antofagasta', '03': 'Atacama', '04': 'Coquimbo', '05': 'Valparaíso', '13': 'Metropolitana', '06': 'O’Higgins', '07': 'Maule', '16': 'Ñuble', '08': 'Biobío', '09': 'Araucanía', '14': 'Los Ríos', '10': 'Los Lagos', '11': 'Aysén', '12': 'Magallanes'}

        data = pd.read_csv(endpoint)
        dates = data.columns[2:]
        
        #['Camas UCI habilitadas','Camas UCI ocupadas COVID-19','Camas UCI ocupadas no COVID-19']
        

        if not type(tstate) == list:
            tstate = [tstate]
            
        counties = [i for i in tstate if len(i)>2 ]
        if counties:
            print('This method doesn\'t support comunas for Chile')
        
        states = [i for i in tstate if len(i)==2]
        statenames = [cuttoname[i] for i in states]

        capacity = []
        occupied_covid = []
        occupied_non_covid = []

        for i in statenames:
            capacity.append(data.loc[((data['Region']==i) & (data['Serie']=='Camas UCI habilitadas'))].iloc[:,2:])
            occupied_covid.append(data.loc[((data['Region']==i) & (data['Serie']=='Camas UCI ocupadas COVID-19'))].iloc[:,2:])
            occupied_non_covid.append(data.loc[((data['Region']==i) & (data['Serie']=='Camas UCI ocupadas no COVID-19'))].iloc[:,2:])
        
            
        capacity =  np.array(capacity).sum(axis=0)[0] 
        occupied_covid =  np.array(occupied_covid).sum(axis=0)[0] 
        occupied_non_covid =  np.array(occupied_non_covid).sum(axis=0)[0] 

        dates = [datetime.strptime(dates[i][:10],'%Y-%m-%d') for i in range(len(dates))]        
        index = np.where(np.array(dates) >= initdate)[0][0] 

        UCI_capacity =list(capacity[index:])
        UCI_use_covid =list(occupied_covid[index:])
        UCI_use_noncovid =list(occupied_non_covid[index:])

        UCI_dates = dates[index:]
        UCI_tr = [(UCI_dates[i]-initdate).days for i in range(len(UCI_dates))]

        # New daily hospitalized
        endpoint2 = 'https://github.com/MinCiencia/Datos-COVID19/raw/master/output/producto91/Ingresos_UCI.csv'
        UCI_d = pd.read_csv(endpoint2)
        UCI_d_dates = UCI_d.columns[1:]
        UCI_d_dates = [datetime.strptime(UCI_d_dates[i][:10],'%Y-%m-%d') for i in range(len(UCI_d_dates))]       
        UCI_d = UCI_d.iloc[0][1:]

        index = np.where(np.array(UCI_d_dates) >= initdate)[0][0] 
        UCI_d = UCI_d[index:]
        UCI_d_dates = UCI_d_dates[index:]
        UCI_d_tr = [(UCI_d_dates[i]-initdate).days for i in range(len(UCI_d_dates))]
       
        if self:
            # Update database
            #print('updating database')
            #self.data = dfappend(self.data,UCI_capacity,UCI_tr,name[0])  
            #self.data = dfappend(self.data,UCI_use_covid,UCI_tr,name[1])  
            #self.data = dfappend(self.data,UCI_use_noncovid,UCI_tr,name[2])              
            self.UCI_capacity = UCI_capacity
            self.UCI_use_covid = UCI_use_covid
            self.UCI_use_noncovid = UCI_use_noncovid                        
            self.UCI_dates = UCI_dates
            self.UCI_tr = UCI_tr
            self.UCI_d = UCI_d
            self.UCI_d_dates = UCI_d_dates
            self.UCI_d_tr = UCI_d_tr
            return
        else:        
            return UCI_capacity,UCI_use_covid,UCI_use_noncovid,UCI_dates,UCI_tr,UCI_d,UCI_d_dates,UCI_d_tr


    # -------------------------------------- #
    #          Deaths (DEIS) MinCiencia      #
    # -------------------------------------- #
    def imp_deaths_deis_mcyt(self=None,tstate = '',initdate = None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto50/DefuncionesDEIS_confirmadosPorComuna.csv',user=None,password=None,name = ['D_confirmed','D_suspected','D_ac_confirmed','D_ac_suspected']):     
        """
            Por ahora solo para data nacional

            Import Accumulated Deaths
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (optional): 
            output: 
                - Br: Real accumulated deaths
                - Br_tr: days from simulation first day
                - Br_dates: data dates
            Usage:
                Br,Br_tr,Br_dates = importAccumulatedDeaths(self=None,tstate = '13',initdate = datetime(2020,5,15))
            
            Compartir la 
        """
        print('Importing Deaths by DEIS')

        if self:
            tstate = self.tstate
            initdate = self.initdate

            
        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")
 
             
        D_r_confirmed = pd.read_csv('https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto50/DefuncionesDEIS_confirmadosPorComuna.csv').dropna().iloc[:,5:].sum(axis=0)
        D_r_suspected = pd.read_csv('https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto50/DefuncionesDEIS_sospechososPorComuna.csv').dropna().iloc[:,5:].sum(axis=0)


        D_r_dates = [datetime.strptime(D_r_confirmed.index[i][:10],'%Y-%m-%d') for i in range(len(D_r_confirmed.index))]  
                 
        index = np.where(np.array(D_r_dates) >= initdate)[0][0] 
        D_r_confirmed = D_r_confirmed[index:]
        D_r_suspected = D_r_suspected[index:]
        D_r_dates = D_r_dates[index:]
        D_r_tr = [(D_r_dates[i]-initdate).days for i in range(len(D_r_dates))]

        B_r_confirmed = D_r_confirmed.cumsum()
        B_r_suspected = D_r_suspected.cumsum()


        if self:
            # Update database
            print('updating database')
            self.data = dfappend(self.data,D_r_confirmed,D_r_tr,name[0])  
            self.data = dfappend(self.data,D_r_suspected,D_r_tr,name[1])  
            self.data = dfappend(self.data,B_r_confirmed,D_r_tr,name[2]) 
            self.data = dfappend(self.data,B_r_suspected,D_r_tr,name[3])             
            
            self.Dr = D_r_confirmed
            self.Dr_suspected = D_r_suspected
            self.Br = B_r_confirmed
            self.Br_suspected = B_r_suspected
            self.Br_dates = D_r_dates
            self.Br_tr = D_r_tr            
            return
        else:        
            return D_r_confirmed,D_r_suspected,B_r_confirmed,B_r_suspected,D_r_dates,D_r_tr


    # ------------------------------------------- #
    #      National Vaccinated Infected Data      #
    # ------------------------------------------- #
    def imp_vacc_infected_mcyt(self=None,initdate = None,endpoint = 'https://github.com/MinCiencia/Datos-COVID19/raw/master/output/producto90/incidencia_en_vacunados.csv',user=None,password=None):     
        """
            Product 90 has data of: 
                * 
            Import Vaccines
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (optional): 
            output: 
                - V_1st: First Dose (2 doses vaccine)
                - V_complete: complete (1 or 2 doses vaccine)
                - V_boost: Boost Dose               
                
            Usage:
                -
        """
        
        print('Importing Vaccines')
        if self:
            initdate = self.initdate
            
        else:
            if not initdate:
                raise Exception("Initial date missing")        

        data = pd.read_csv(endpoint) 
        #data['date'] = [datetime(2021,1,3) + timedelta(days=7*(i-1)) for i in data['semana_epidemiologica']]
        data['date'] = [datetime(2021,1,4) + timedelta(days=7*(i)) for i in data['semana_epidemiologica'].index]
        
        initweek = int((initdate - datetime(2021,1,3)).days/7)    

        Iv_df = data[['dos_dosis_comp_casos','dosis_unica_comp_casos','dosis_ref_comp_casos']]
        Inv_df = data['sin_vac_casos']
        
        Iv_s = data['dos_dosis_comp_casos'] + data['dosis_unica_comp_casos'] + data['dosis_ref_comp_casos']        
        

        Iv_d = []
        Inv_d = []
        
        
        for i in Iv_s:
            for j in range(7):
                Iv_d.append(round(i/7))

        for i in Inv_df:
            for j in range(7):
                Inv_d.append(round(i/7))

        Iv_dates = [datetime(2021,1,4)+timedelta(days=i) for i in range(len(Iv_d))]     

        Iv_d = pd.DataFrame({'Iv_d':Iv_d},index=Iv_dates).rolling(window=7,win_type='gaussian', min_periods=1, center=True).mean(std=2).round()
        Inv_d = pd.DataFrame({'Inv_d':Inv_d},index=Iv_dates).rolling(window=7,win_type='gaussian', min_periods=1, center=True).mean(std=2).round()

        Iv_ac = Iv_d.cumsum()
        Inv_ac = Inv_d.cumsum()
        
        Iv_d = Iv_d.loc[np.array(Iv_dates)>initdate]
        Iv_ac = np.array(Iv_ac.loc[np.array(Iv_dates)>initdate]['Iv_d'])
        Inv_d = Inv_d.loc[np.array(Iv_dates)>initdate]
        Inv_ac = np.array(Inv_ac.loc[np.array(Iv_dates)>initdate]['Inv_d'])

        Iv_dates = Iv_dates[np.where(np.array(Iv_dates) > initdate)[0][0]:]
    
        # Active vaccinated infected
        Iv0 = data['dos_dosis_comp_casos'][initweek-1:initweek+1].sum() + data['dosis_unica_comp_casos'][initweek-1:initweek+1].sum() + data['dosis_ref_comp_casos'][initweek-1:initweek+1].sum() 
        Iv_d = np.array(Iv_d['Iv_d'])
        Inv_d = np.array(Inv_d['Inv_d'])
        
        if self:
            self.Iv_ac = Iv_ac
            self.Iv_d = Iv_d
            self.Iv0 = Iv0
            self.Inv_ac = Inv_ac
            self.Inv_d = Inv_d            
            self.Iv_dates = Iv_dates
            self.Iv_df = Iv_df
            return
        else:        
            return Iv_ac,Iv_d,Iv0, Inv_d,Inv_ac,Iv_dates


    # ------------------------------------------- #
    #      National Vaccinated Data      #
    # ------------------------------------------- #
    def imp_vaccinated_mcyt(self=None,tstate = '',initdate = None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto76/vacunacion.csv',user=None,password=None):     
        """
            Product 76 has data of: 
                * 
            Import Vaccines
            input: 
                - tstate: [string or string list] CUT por comuna o región
                - initdate: datetime object with the initial date
                - endpoint (optional): 
            output: 
                - V_1st: First Dose (2 doses vaccine)
                - V_complete: complete (1 or 2 doses vaccine)
                - V_boost: Boost Dose               
                
            Usage:
                -
        """
        tstate = ['15','01','02','03','04','05','13','06','07','16','08','09','14','10','11','12'] 

        print('Importing Vaccines')
        if self:
            tstate = self.tstate
            initdate = self.initdate
            population = self.population

        else:
            if not tstate:
                raise Exception("State code missing")
            if not initdate:
                raise Exception("Initial date missing")
            population = 19120000

        data = pd.read_csv(endpoint)
        
        # Vacunas acumuladas
        v1_ac = data.loc[(data['Region']=='Total')&(data['Dosis']=='Primera')].iloc[:,2:].iloc[0]
        v2_ac = data.loc[(data['Region']=='Total')&(data['Dosis']=='Segunda')].iloc[:,2:].iloc[0]
        v3_ac = data.loc[(data['Region']=='Total')&(data['Dosis']=='Refuerzo')].iloc[:,2:].iloc[0]

        # Vacunas diarias:
        v1_d = v1_ac.diff()
        v2_d = v2_ac.diff()
        v3_d = v3_ac.diff()
        
        v1_d[0] = v1_ac[0]
        v2_d[0] = v1_ac[0]
        v3_d[0] = v1_ac[0]

        # Personas activas en estado de n dosis
        v3 = v3_ac
        v2 = v2_ac-v3_ac
        v1 = v1_ac-v2_ac
        v0 = population - v1_ac

        v_dates = [datetime.strptime(data.columns[2:][i],'%Y-%m-%d') for i in range(len(data.columns[2:]))]

        # Select dates: 
        index = np.where(np.array(v_dates) >= initdate)[0][0]

        v1_d = v1_d[index:]
        v2_d = v2_d[index:]
        v3_d = v3_d[index:]

        v1_ac = v1_ac[index:]
        v2_ac = v2_ac[index:]
        v3_ac = v3_ac[index:]

        v1 = v1[index:]
        v2 = v2[index:]
        v3 = v3[index:]

        v1.name = '1 Dose Active'
        v2.name = '2 Doses Active'
        v3.name = '3 Doses Active'

        v1_d.name = '1 Dose application'
        v2_d.name = '2 Doses application'
        v3_d.name = '3 Doses application'

        v1_ac.name = '1 Dose accumulated'
        v2_ac.name = '2 Doses accumulated'
        v3_ac.name = '3 Doses accumulated'        
        
        if self:
            self.v1_d = v1_d
            self.v2_d = v2_d
            self.v3_d = v3_d

            self.v1_ac = v1_ac
            self.v2_ac = v2_ac
            self.v3_ac = v3_ac

            self.v1 = v1
            self.v2 = v2
            self.v3 = v3

            self.v_dates = v_dates
            return
        else:
            return v1,v2,v3,v1_d,v2_d,v3_d,v1_ac,v2_ac,v3_ac, v_dates

