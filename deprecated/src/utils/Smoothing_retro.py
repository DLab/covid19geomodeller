# -*- coding: utf-8 -*-
import pandas as pd
import datetime
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

'''
Ejemplo de suavizado de curva con retropropagacion de datos excesivos

Author: Felipe Castillo
fcastillo@dlab.cl

'''


def prepare_cases_retro(cases, cutoff=0):
    '''
    input: interpolated cases
    output: original interpolated cases, smoothed serie using a gaussian moving window

    se remueven los ~30K entre las fechas en que se interpolaron esos datos y luego se retropropagan proporcionalmente
    '''

    new_cases = cases.diff() # pasa de acumulados a diarios
    new_cases[datetime.date(2020,6,16):datetime.date(2020,6,19)]=new_cases[datetime.date(2020,6,16):datetime.date(2020,6,19)]-31412/3 #quitando los ~30K de donde fueron interpolados
    new_cases[:datetime.date(2020,6,19)] += np.round(new_cases[:datetime.date(2020,6,19)]*31412/np.sum(new_cases[:datetime.date(2020,6,19)])) #Retropropagando los ~30K
    smoothed = new_cases.rolling(7, win_type='gaussian', min_periods=1, center=True).mean(std=2).round()
    idx_start = np.searchsorted(smoothed, cutoff)

    smoothed = smoothed.iloc[idx_start:]
    original = new_cases.loc[smoothed.index]
    original = new_cases

    return original, smoothed




### Cargando datos y dropeando columnas innecesarias
endpoint = 'https://raw.githubusercontent.com/MinCiencia/Datos-COVID19/master/output/producto1/Covid-19.csv'
comunas = pd.read_csv(endpoint)
comunas = comunas.fillna(0)

comunas = comunas.drop(['Tasa', 'Poblacion'], axis=1)
comunas.rename(columns = {'Codigo region':'id_region', 'Codigo comuna':'id_comuna'}, inplace=True)

comunas.index = comunas['Comuna']
dates = comunas.drop(columns=['Region', 'id_region', 'Comuna', 'id_comuna']).columns

regiones = comunas.pivot_table(index='Region', values=dates, aggfunc=sum)
regiones = regiones.sum(axis=0)
regiones.index = pd.to_datetime(regiones.index)

### Interpolación

inicio = datetime.date(2020,4,1)
fin = datetime.date(2020,regiones.index[-1].month,regiones.index[-1].day)
dias = [inicio + datetime.timedelta(days=d) for d in range((fin - inicio).days + 1)]

### Data frame para almacenar la interpolacion
data = pd.DataFrame(index=dias)
data.loc[regiones.index[1:],'Nacional'] = regiones.loc[datetime.datetime(2020,4,1):].values



### Interpolación
t = np.arange(0,len(dias),1)
vec_ = np.array(list(data.loc[datetime.date(2020,3,30):datetime.date(2020,data.index[-1].month,data.index[-1].day +1)].values))
interp = interp1d(t[~np.isnan(vec_)[:,0]], vec_[~np.isnan(vec_)])
interp_lineal=interp(t)

data['Nacional'] = interp_lineal

# Retropropagando los 30K
### Smoothing

new, smoothed = prepare_cases_retro(data['Nacional'])

### PLOT
plt.plot(smoothed)
plt.show()
