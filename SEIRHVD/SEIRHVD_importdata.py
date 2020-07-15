#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import requests
import pandas as pd
from datetime import datetime
import numpy as np

"""   
   SEIRHDV Import data
"""

class SEIRHVD_importdata():
    # ------------------------------- #
    #        Importar Data Real       #
    # ------------------------------- #

    def importinfectadosactivos(self):
        # ---------------------- # 
        #   Infectados Activos   #
        # ---------------------- #
        cutlist = []
        cutlistpath = "../Data/cutlist.csv"
        cutlist = pd.read_csv(cutlistpath, header = None,dtype=str)

        actives = []
        mydict = None
        for index, row in cutlist.iterrows():    
            state = str(row[0])[0:2]
            comuna = str(row[0])
            if self.tstate == state:
                endpoint = "http://192.168.2.223:5006/getActiveNewCasesByComuna?comuna="+comuna
                r = requests.get(endpoint) 
                mydict = r.json()
                actives.append(mydict['actives'])
                #data=pd.DataFrame(mydict)
        self.Ir = (np.array(actives)).sum(axis=0)
        self.Ir_dates = [datetime.strptime(mydict['dates'][i],'%Y-%m-%d') for i in range(len(mydict['dates']))]

        index = np.where(np.array(self.Ir_dates) >= self.initdate)[0][0]     
        self.Ir=self.Ir[index:]
        self.Ir_dates=self.Ir_dates[index:]
        self.tr = [(self.Ir_dates[i]-self.initdate).days for i in range(len(self.Ir))]
        print('Infectados Activos')
        return



    # -------------------------------- #
    #    Datos Infectados acumulados   #
    # -------------------------------- #
    def importinfectadosacumulados(self,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv' ):     
        aux = pd.read_csv(endpoint)        
        self.I_ac_r = aux.loc[aux['Codigo region']==int(self.tstate)].iloc[:,5:-1].sum()
        
        self.I_ac_r_dates = [datetime.strptime(self.I_ac_r.index[i],'%Y-%m-%d') for i in range(len(self.I_ac_r))]
        index = np.where(np.array(self.I_ac_r_dates) >= self.initdate)[0][0] 
        self.I_ac_r = self.I_ac_r[index:]
        self.I_ac_r_dates = self.I_ac_r_dates[index:]
        self.I_ac_r_tr = [(self.I_ac_r_dates[i]-self.initdate).days for i in range(len(self.I_ac_r))]
        print('Infectados Acumulados')
        return



    # -------------------------------- #
    #    Datos Infectados diarios   #
    # -------------------------------- #
    def importinfectadosdiarios(self,endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv' ):     
        aux = pd.read_csv(endpoint)        
        self.I_d_r = aux.loc[aux['Codigo region']==int(self.tstate)].iloc[:,5:-1].sum().diff()
        
        self.I_d_r_dates = [datetime.strptime(self.I_d_r.index[i],'%Y-%m-%d') for i in range(len(self.I_d_r))]
        index = np.where(np.array(self.I_d_r_dates) >= self.initdate)[0][0] 
        self.I_d_r = self.I_d_r[index:]
        self.I_d_r_dates = self.I_d_r_dates[index:]
        self.I_d_r_tr = [(self.I_d_r_dates[i]-self.initdate).days for i in range(len(self.I_d_r))]
        print('Infectados diarios')
        return

    def importsochimi(self,endpoint = "http://192.168.2.223:5006/getBedsAndVentilationByState?state="):
        # ------------------ #
        #    Datos Sochimi   #
        # ------------------ #
        endpoint = endpoint+self.tstate
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
        cut =  ['15','01','02','03','04','05','13','06','07','16','08','09','14','10','11','12','00']
        index = cut.index(self.tstate)
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
        index = cut.index(self.tstate)
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
        self.ED_RM_df = excess_dead.loc[excess_dead['Codigo region']==int(self.tstate)]      
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

    # --------------------------- #
    #    Importar toda la data    #
    # --------------------------- #

    def importdata(self):
        print('Importando Datos')
        self.importfallecidosacumulados()
        self.importfallecidosexcesivos()
        self.importinfectadosactivos()
        self.importsochimi()
        self.importpcrpop()
        self.importinfectadosdiarios()
        self.importinfectadosacumulados()
        print('Done')
