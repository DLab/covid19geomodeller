
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from datetime import datetime
from datetime import timedelta
import pandas as pd
import numpy as np
"""
# -------------------------- #
#    Generación de Tablas    #
# -------------------------- #


# Variables: 
# I_cum: Infectados acumulados
# I_act: Infectados activos
# D: Muertos acumulados
# D_d: Muertos Diarios
# L: Letalidad
# H: Uso de Camas hospitalarias
# V: Uso de VMI
# H_tot: Necesidad total de camas (incluidas las que se necesitan por sobre la capacidad)
# V_tot: Necesidad total de VMI (incluidos las que se necesitan por sobre la capacidad)


"""
class SEIRHVD_tables():    
    # ------------------------------------- #
    #    Results data table construction    #
    # ------------------------------------- #


    #from datetime import timedelta
    def tabladedatos(self,inicio = datetime(2020,5,15), fin = datetime(2020,6,30),variables =['I_cum','I_act','D','L'], path=''):

        # Time
        tr_i = (inicio-self.initdate).days
        tr_f = (fin-self.initdate).days
        days = (fin-inicio).days
        index = [(inicio+timedelta(days=i)).strftime("%d/%m/%Y") for i in range(days+1)]
        idx = [np.searchsorted(self.t[i],range(tr_i,tr_f+1)) for i in range(self.numescenarios)]

        #data
        data = []
        # --------------------------- #
        #    Fallecidos acumulados    #
        # --------------------------- #
        if 'D' in variables:
            # Interpolacion para datos faltantes
            Bdata = dict()
            namelist = ['Muertos-60','Muertos-65','Muertos-70']
            for i in range(len(namelist)):        
                B_hoy = [round(self.B[i][idx[i][0]])]
                for j in range(1,len(idx[i])-1):
                    if idx[i][j-1]== idx[i][j]:
                        B_hoy.extend([round((self.B[i][idx[i][j-1]]+self.B[i][idx[i][j+1]])/2)])
                    else:        
                        B_hoy.extend([round(self.B[i][idx[i][j]])])
                B_hoy.extend([round(self.B[i][idx[i][-1]])])
                Bdata[namelist[i]]=B_hoy

            Bdata = pd.DataFrame(Bdata)
            data.append(Bdata)


        # Infectados Acumulados
        if 'I_cum' in variables:
            Iacdata = dict()
            namelist = ['Infectados-60','Infectados-65','Infectados-70']
            for i in range(self.numescenarios):        
                Iac_hoy = [round(self.Iac[i][idx[i][0]])]
                for j in range(1,len(idx[i])-1):
                    if idx[i][j-1]== idx[i][j]:
                        Iac_hoy.extend([round((self.Iac[i][idx[i][j-1]]+self.Iac[i][idx[i][j+1]])/2)])
                    else:        
                        Iac_hoy.extend([round(self.Iac[i][idx[i][j]])])
                Iac_hoy.extend([round(self.Iac[i][idx[i][-1]])])
                Iacdata[namelist[i]]=Iac_hoy
            Iacdata = pd.DataFrame(Iacdata)
            data.append(Iacdata)


        # Letalidad
        #let = [100*B[i]/Iac[i] for i in range(self.numescenarios)
        #Iacdata = dict()
        #namelist = ['Infectados-60','Infectados-65','Infectados-70']
        #for i in range(numescenarios):        
        #    Iac_hoy = [round(Iac[i][idx[i][0]])]
        #    for j in range(1,len(idx[i])-1):
        #        if idx[i][j-1]== idx[i][j]:
        #            Iac_hoy.extend([round((Iac[i][idx[i][j-1]]+Iac[i][idx[i][j+1]])/2)])
        #        else:        
        #            Iac_hoy.extend([round(Iac[i][idx[i][j]])])
        #    Iac_hoy.extend([round(Iac[i][idx[i][-1]])])
        #    Iacdata[namelist[i]]=Iac_hoy
        #
        #Iacdata = pd.DataFrame(Iacdata)




        # ------------------ #
        #    Uso de Camas    #
        # ------------------ #

        #H_bed_hoy = [H_bed[i][idx[i][0:days]] for i in range(len(input))]
        if 'H' in variables:
            UsoCamas = dict()
            namelist = ['UsoCamas-60','UsoCamas-65','UsoCamas-70']
            for i in range(3):        
                UsoCamas_hoy = [round(self.H_bed[i][idx[i][0]])]
                for j in range(1,len(idx[i])-1):
                    if idx[i][j-1]== idx[i][j]:
                        UsoCamas_hoy.extend([round((self.H_bed[i][idx[i][j-1]]+self.H_bed[i][idx[i][j+1]])/2)])
                    else:        
                        UsoCamas_hoy.extend([round(self.H_bed[i][idx[i][j]])])
                UsoCamas_hoy.extend([round(self.H_bed[i][idx[i][-1]])])
                UsoCamas[namelist[i-3]]=UsoCamas_hoy

            Hbed = pd.DataFrame(UsoCamas)
            data.append(Hbed)

        # ------------------------- #
        #    Uso de Ventiladores    #
        # ------------------------- #
        #H_vent_hoy = [H_vent[i][idx[i][0:days]] for i in range(len(input))]
        if 'V' in variables:
            namelist = ['UsoVMI-60','UsoVMI-65','UsoVMI-70']
            UsoVMI = dict()
            for i in range(3):        
                UsoVMI_hoy = [round(self.H_vent[i][idx[i][0]])]
                for j in range(1,len(idx[i])-1):
                    if idx[i][j-1]== idx[i][j]:
                        UsoVMI_hoy.extend([round((self.H_vent[i][idx[i][j-1]]+self.H_vent[i][idx[i][j+1]])/2)])
                    else:        
                        UsoVMI_hoy.extend([round(self.H_vent[i][idx[i][j]])])
                UsoVMI_hoy.extend([round(self.H_vent[i][idx[i][-1]])])
                UsoVMI[namelist[i-3]]=UsoVMI_hoy

            Hvent = pd.DataFrame(UsoVMI)
            data.append(Hvent)


        # ---------------------------------- #
        #    Camas adicionales Requeridas    #
        # ---------------------------------- #
        if 'H_ad' in variables:
            CH_d = dict()
            namelist = ['CamaAdicional-60','CamaAdicional-65','CamaAdicional-70']
            for i in range(3):        
                CH_hoy = [round(self.CH[i][idx[i][0]])]
                for j in range(1,len(idx[i])-1):
                    if idx[i][j-1]== idx[i][j]:
                        CH_hoy.extend([round((self.CH[i][idx[i][j-1]]+self.CH[i][idx[i][j+1]])/2)])
                    else:        
                        CH_hoy.extend([round(self.CH[i][idx[i][j]])])
                CH_hoy.extend([round(self.CH[i][idx[i][-1]])])
                CH_d[namelist[i-3]]=CH_hoy

            CH_d = pd.DataFrame(CH_d)
            data.append(CH_d)

        if 'V_ad' in variables:
            CV_d = dict()
            namelist = ['VMIAdicional-60','VMIAdicional-65','VMIAdicional-70']
            for i in range(3):        
                CV_hoy = [round(self.CV[i][idx[i][0]])]
                for j in range(1,len(idx[i])-1):
                    if idx[i][j-1]== idx[i][j]:
                        CV_hoy.extend([round((self.CV[i][idx[i][j-1]]+self.CV[i][idx[i][j+1]])/2)])
                    else:        
                        CV_hoy.extend([round(self.CV[i][idx[i][j]])])
                CV_hoy.extend([round(self.CV[i][idx[i][-1]])])
                CV_d[namelist[i-3]]=CV_hoy

            CV_d = pd.DataFrame(CV_d)
            data.append(CV_d)


        # ------------------------------ #
        #    Necesidad total de Camas    #
        # ------------------------------ #
        if False:
            namelistUsoC = ['UsoCamas-60','UsoCamas-65','UsoCamas-70']
            namelistUsoV = ['UsoVMI-60','UsoVMI-65','UsoVMI-70']

            namelistCH = ['CamaAdicional-60','CamaAdicional-65','CamaAdicional-70']
            namelistCV = ['VMIAdicional-60','VMIAdicional-65','VMIAdicional-70']

            namelistcamas = ['NecesidadTotalCamas-60','NecesidadTotalCamas-65','NecesidadTotalCamas-70']
            namelistvmi = ['NecesidadTotalVMI-60','NecesidadTotalVMI-65','NecesidadTotalVMI-70']

            totbed_d = pd.DataFrame()

            totvmi_d = pd.DataFrame()

            for i in range(len(namelistCH)):
                totbed_d[namelistcamas[i]] = CH_d[namelistCH[i]] + Hbed[namelistUsoC[i]]
                totvmi_d[namelistvmi[i]] = CV_d[namelistCV[i]] + Hvent[namelistUsoV[i]]
            data.append(totbed_d)
            data.append(totvmi_d )


        # ------------------------- #
        #     Create Data Frame     #
        # ------------------------- #
        index = pd.DataFrame(dict(dates=index))
        data = pd.concat(data, axis=1, sort=False)
        data = pd.concat([index,data], axis=1, sort=False) 
        data = data.set_index('dates')
        if path:
            data.to_excel(path)
            #
        return(data)
