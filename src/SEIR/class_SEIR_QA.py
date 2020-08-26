

"""

 SEIR tests debug y QA

"""

import sys
from pathlib import Path
sys.path.insert(1, '../SEIR/')
sys.path.insert(1, 'SEIR/')


from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing

import numpy as np
import pandas as pd
import requests

from class_SEIR import SEIR
from Quarantine import Quarantine


tsim = 1000
max_mob = 0.85
rem_mob = 0.4 
qp = 0
iqt = 0
alpha = Quarantine.alphafunct(rem_mob,iqt = iqt)

beta = 0.117
mu = 1.5
k = 40
I_ac = 100
I = 100
population = 1000000

model = SEIR(tsim,alpha,beta,mu,k=0,I=I,I_ac=I_ac,population=population)

# Integrate
model.integr_sci(0,tsim,0.1)
model.integr(0,tsim,0.1)

"""

   QA

"""
# -------------------------- #
#        Plot function       #
# -------------------------- #
def plot(title = '',xlabel='',ylabel='',legend=True):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if legend:
        plt.legend(loc=0)
    plt.show()
    

Delta_population = model.S+model.E+model.I+model.R - population
max(Delta_population)


# Plots
plt.plot(model.t,model.S,label='S')
plt.plot(model.t,model.E,label='E')
plt.plot(model.t,model.I,label='I')
plt.plot(model.t,model.R,label='R')
plot(title = 'SEIR - QA')




"""
# ------------------------------ #

          Meta Simulation

# ------------------------------ #
"""


# ------------------------------ #
#           Import data          #
# ------------------------------ #

tstate = '13'
initdate = datetime(2020,5,15)


# Infectados Acumulados y Diarios
endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv'

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


# Infectados Activos

cutlist = []
cutlistpath = "Data/cutlist.csv"
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

# Population:

endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv'
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



# ------------------------------ #
#           Simulation           #
# ------------------------------ #

from SEIR_parallel import seirMetaAnalysis


s1 = Quarantine.alphafunct(0.5)
s2 = Quarantine.alphafunct(0.6)
s3 = Quarantine.alphafunct(0.7)
quarantines = [s1,s2,s3] # Cambar a quarantines


k = [0,5,10,15,20]

meta = seirMetaAnalysis()

# Simulate for k
sims = meta.simulate_k2(tsim,quarantines,beta,mu,k=k, I=Ir[0],I_ac=I_ac_r[0],I_d=I_d_r[0],population=population )

# ------------------------------ #
#            Analysis            #
# ------------------------------ #
realdata = True
xlim = 100

# Actives
ylim = max([max(sims[i][j].I[:xlim]) for j in range(np.shape(sims)[1]) for i in range(np.shape(sims)[0]) ])
fig, axs = plt.subplots(len(k), len(quarantines))

for i in range(len(k)):
    for j in range(len(quarantines)):
        axs[i, j].plot(sims[i][j].t,sims[i][j].I,label="Infectados")
        axs[i, j].set_title("K: "+str(k[i])+" | Alpha: "+str([0.5,0.6,0.7][j]))
        if realdata == True:
            axs[i, j].scatter(tr,Ir,label='Infectados Activos reales')        
        axs[i, j].set_ylim([0,ylim*1.05])
        axs[i, j].set_xlim([0,xlim])             

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Infectados Activos')
fig.show()


# Acumulados
ylim = max([max(sims[i][j].I_ac[:xlim]) for j in range(np.shape(sims)[1]) for i in range(np.shape(sims)[0]) ])
fig, axs = plt.subplots(len(k), len(quarantines))

for i in range(len(k)):
    for j in range(len(quarantines)):
        axs[i, j].plot(sims[i][j].t,sims[i][j].I_ac,label="Infectados")
        axs[i, j].set_title("K: "+str(k[i])+" | Alpha: "+str([0.5,0.6,0.7][j]))
        if realdata == True:
            axs[i, j].scatter(I_ac_r_tr,I_ac_r,label='Infectados Acumulados reales')
        axs[i, j].set_ylim([0,ylim*1.05])
        axs[i, j].set_xlim([0,xlim])             

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Infectados Acumulados')
fig.show()


# Diarios
ylim = max([max(sims[i][j].I_d[:xlim]) for j in range(np.shape(sims)[1]) for i in range(np.shape(sims)[0]) ])
fig, axs = plt.subplots(len(k), len(quarantines))

for i in range(len(k)):
    for j in range(len(quarantines)):
        axs[i, j].plot(sims[i][j].t,sims[i][j].I_d,label="Infectados")
        axs[i, j].set_title("K: "+str(k[i])+" | Alpha: "+str([0.5,0.6,0.7][j]))
        if realdata == True:
            axs[i, j].scatter(I_d_r_tr,I_d_r,label='Infectados diarios reales')        
        axs[i, j].set_ylim([0,ylim*1.05])
        axs[i, j].set_xlim([0,xlim])             

lines, labels = fig.axes[-1].get_legend_handles_labels()  
fig.legend(lines, labels,loc = 'best')
fig.suptitle('Infectados Diarios')
fig.show()

