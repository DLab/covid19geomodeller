import sys
from pathlib import Path
sys.path.insert(1, '../src/SEIRHVD/')
sys.path.insert(1, 'src/SEIRHVD/')
sys.path.insert(1, '../src/utils/')
sys.path.insert(1, 'src/utils/')
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd

from SEIRHVD_local import SEIRHVD_local
from Quarantine import Quarantine


tstate = '05'
# Fecha Inicial
initdate = datetime(2020,5,15)
# Parametros del modelo
beta = 0.2 # Tasa de contagio
mu = 0.6 # Razon E0/I0
ScaleFactor = 1 # Factor de Escala: Numero de infectados por sobre los reportados
SeroPrevFactor = 1 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos
tsim = 500 # Tiempo de simulacion
# Creación del objeto de simulación
simulation = SEIRHVD_local(beta = beta,mu = mu,ScaleFactor=ScaleFactor,SeroPrevFactor=SeroPrevFactor,expinfection=expinfection,initdate = initdate, tsim = tsim,tstate=tstate)
quarantines = [[500.0, 0.85, 0.7, 0.0, 0.0, 500.0, 0.0]]
simulation.inputarray = np.array(quarantines)
simulation.addquarantine()
# Valores iniciales
I_act0 = 100
dead = 0
population = 10000
H0 = 1
V0 = 1
# Capacidades hospitalarias
Htot = 1000
Vtot = 1000
simulation.initialvalues(I_act0,dead,population,H0,V0,Htot,Vtot,R=0,D=0,H_cr = 0)
# Simular
simulation.simulate()
simulation.plotinfectadosactivos()




# Import Quarantine data}
def importQuarantines(tstate):
    """
    Import historical quarantines
    tstate could be a region, a comuna, or a list with a combination of both.
    """
    endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto29/Cuarentenas-Totales.csv'
    aux = pd.read_csv(endpoint)        
    aux['Codigo CUT Comuna'] = aux['Código CUT Comuna'].map(lambda i: str(i).zfill(5))
    aux = aux.loc[[tstate in x[:len(tstate)] for x in aux['Codigo CUT Comuna']]]
    
    if type(tstate) == list:
        cuarentenas = pd.DataFrame()
        for i in tstate:
            cuarentenas = cuarentenas.append(aux.loc[[i in x[:len(i)] for x in aux['Codigo CUT Comuna']]])
    else:                        
        cuarentenas = aux.loc[[tstate in x[:len(tstate)] for x in aux['Codigo CUT Comuna']]]

    return cuarentenas


# Valpo
tstate = '05'
# Bio Bio
tstate = '08' 
# Araucania
tstate = '09' 

# Import Data
initdate = datetime(2020,5,15)
sochimi,sochimi_dates ,sochimi_tr,Hr, Hr_tot, Vr ,Vr_tot = a.importsochimi(tstate=tstate,initdate=initdate)   
I_d_r_smooth, I_d_r_tr, I_d_r_dates = a.importDailyInfected(tstate=tstate,initdate = initdate) 

# Plot Data

# Infected
plt.plot(I_d_r_dates,I_d_r_smooth)
plt.title('Daily Infected - Region '+tstate)
plt.show()

# Bed use
plt.plot(sochimi_dates,Hr,label='Bed Use')
plt.plot(sochimi_dates,Hr_tot,label='Bed Capacity')
plt.plot(sochimi_dates,Vr,label='VMI use')
plt.plot(sochimi_dates,Vr_tot,label='VMI Capacity')
plt.title('SOCHIMI - Region '+tstate)
plt.legend(loc=0)
plt.show()
