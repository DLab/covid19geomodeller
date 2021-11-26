#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import requests
import pandas as pd
from datetime import datetime
from datetime import timedelta
import numpy as np

"""   
   SEIRHDV Import data
"""

class SEIR_importdata():
    # ------------------------------- #
    #        Importar Data Real       #
    # ------------------------------- #

    # -------------------------------- #
    #            Poblacion             #
    # -------------------------------- #
    def importPopulation(self=None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv',tstate = ''):     
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

    def importinfectadosactivos(self=None,tstate = '',initdate=None):
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
        cutlist = []
        cutlistpath = "../Data/cutlist.csv"
        cutlist = pd.read_csv(cutlistpath, header = None,dtype=str)

        actives = []
        mydict = None
        if type(tstate) == list:
            for i in tstate:
                if len(i)==2:
                    for index, row in cutlist.iterrows():    
                        state = str(row[0])[0:2]
                        comuna = str(row[0])
                        if i == state:
                            endpoint = "http://192.168.2.223:5006/getActiveNewCasesByComuna?comuna="+comuna
                            r = requests.get(endpoint) 
                            mydict = r.json()
                            actives.append(mydict['actives'])
                            #data=pd.DataFrame(mydict)
                    #Ir = (np.array(actives)).sum(axis=0)
                elif len(i)>2:
                    endpoint = "http://192.168.2.223:5006/getActiveNewCasesByComuna?comuna="+i
                    r = requests.get(endpoint) 
                    mydict = r.json()
                    actives.append(mydict['actives'])
                    #Ir = np.array(mydict['actives'])
                Ir = (np.array(actives)).sum(axis=0)
        else:
            if len(tstate)==2:
                for index, row in cutlist.iterrows():    
                    state = str(row[0])[0:2]
                    comuna = str(row[0])
                    if tstate == state:
                        endpoint = "http://192.168.2.223:5006/getActiveNewCasesByComuna?comuna="+comuna
                        r = requests.get(endpoint) 
                        mydict = r.json()
                        actives.append(mydict['actives'])
                        #data=pd.DataFrame(mydict)
                Ir = (np.array(actives)).sum(axis=0)
            elif len(tstate)>2:
                endpoint = "http://192.168.2.223:5006/getActiveNewCasesByComuna?comuna="+tstate
                r = requests.get(endpoint) 
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



    # -------------------------------- #
    #    Datos Infectados acumulados   #
    # -------------------------------- #
    def importinfectadosacumulados(self=None,tstate = '',initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv'):     
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



    # -------------------------------- #
    #    Datos Infectados diarios      #
    # -------------------------------- #
    # Falta interpolar
    def importinfectadosdiarios(self=None,tstate = '',initdate = None,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv' ):     

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

        # ------------------ #
        #    Datos Sochimi   #
        # ------------------ #
    def importsochimi(self=None,tstate = '', initdate = None, endpoint = "http://192.168.2.223:5006/getBedsAndVentilationByState?state="):

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
        self.sochimi=pd.DataFrame(mydict)
        sochimi = self.sochimi
        self.Hr = sochimi['camas_ocupadas']
        self.Vr =  sochimi['vmi_ocupados']
        self.Vr_tot =  sochimi['vmi_totales']
        self.Hr_tot =  sochimi['camas_totales']
        self.sochimi_dates = [datetime.strptime(sochimi['dates'][i][:10],'%Y-%m-%d') for i in range(len(sochimi))]

        index = np.where(np.array(self.sochimi_dates) >= self.initdate)[0][0] 
        self.Hr=list(self.Hr[index:])
        self.Vr=list(self.Vr[index:])
        self.Hr_tot=list(self.Hr_tot[index:])
        self.Vr_tot=(list(self.Vr_tot[index:]))
        self.sochimi_dates = self.sochimi_dates[index:]
        self.sochimi_tr = [(self.sochimi_dates[i]-self.initdate).days for i in range(len(self.Hr))]
        print('Sochimi')
        return(sochimi)

    # -------------------------------- #
    #    Datos Fallecidos acumulados   #
    # -------------------------------- #
    def importfallecidosacumulados(self,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto14/FallecidosCumulativo.csv' ):     
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
        self.Br = pd.read_csv(endpoint).iloc[index][1:] 
        self.Br_dates = [datetime.strptime(self.Br.index[i],'%Y-%m-%d') for i in range(len(self.Br))]
        index = np.where(np.array(self.Br_dates) >= self.initdate)[0][0] 
        self.Br = self.Br[index:]
        self.Br_dates = self.Br_dates[index:]
        self.Br_tr = [(self.Br_dates[i]-self.initdate).days for i in range(len(self.Br))]
        print('Fallecidos Acumulados')
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


    # ---------------------------------------- #
    #    Datos Infectados activos Minciencia   #
    # ---------------------------------------- #
    def importinfectadosactivosminciencia(self=None,tstate = '', initdate = None, endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto19/CasosActivosPorComuna.csv' ):     
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
            return I_minciencia_r, I_minciencia_r_dates, I_minciencia_r_tr


    # ---------------------------------------- #
    #       Datos Subreporte de Infectados     #
    # ---------------------------------------- #
    def importSubreporte(self = None,tstate = '', initdate = None):
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
        

    # --------------------------- #
    #    Importar toda la data    #
    # --------------------------- #

    def importdata(self):
        print('Importando Datos')
        self.importfallecidosacumulados()
        #self.importfallecidosexcesivos()        
        self.importinfectadosactivosminciencia()
        self.importsochimi()
        #self.importpcrpop()
        self.importPopulation()
        self.importinfectadosactivos()
        self.importinfectadosdiarios()
        self.importinfectadosacumulados()
        print('Done')