import numpy as np
import pandas as pd
import toml
from datetime import datetime
from datetime import timedelta


import os
import sys
path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(1, path)

import data.cv19data as cv19data
import utils.cv19timeutils as cv19timeutils
import utils.cv19functions as cv19functions
import utils.cv19files as cv19files
from copy import deepcopy

#from importlib import import_module

"""
Todo:
    [x] Resolver la interaccion del cv19functions con las listas    
    [x] Decidir como decidir qué modelo usar y luego como importarlo
    [x] Encontrar los parámetros "iterables" (listas) en el config y el kwargs
    [x] Crear función recursiva que cree for anidados que recorran los elementos del los parámetros iterables
    [x] Construir arreglo de diccionarios de configuración
    [ ] Instanciar objetos del modelo correspondiente a partir del arreglo de configuración
    [ ] Crear función que simule las funciones en todos los arreglos. Fundamental paralelizar! 

    [ ] Construir una función resumen que imprima las características principales
        * tipo de modelo
        * variables a iterar
        * Condiciones iniciales reales 
"""

class cv19sim():
    def __init__(self,config,model='SEIR', inputdata = None, **kwargs):

        config = deepcopy(config)            
        if model == 'SEIR':
            from models.SEIR import SEIR
            model = SEIR
        elif model == 'SEIRHVD':
            from models.SEIRHVD import SEIRHVD
            model = SEIRHVD

        # Leer el archivo de configuracion
        if not type(config) == dict:
            config = toml.load(config)


        sims = []
        iterables = {}
        # create auxiliar object         
        aux = cv19files.unwrapconfig(config,**kwargs)        
        
        # Find iterable parameters:
        for key,value in aux.items():
            if type(value)==list:
                print(key+':'+str(value))
                iterables.update({key:value})

        #print('There are '+ str(len(iterables))+' iterable parameters')
               

        if iterables:
            # Pop iterables from kwargs
            for key in iterables:
                if key in kwargs:
                    kwargs.pop(key)
            expandediterables = iterate(iterables)
            buildsim = simapply(config,model,inputdata,**kwargs)            
            self.sims = buildsim(expandediterables)
        else:
            self.sims = [model(config,inputdata,**kwargs)]
            print('Simulating over 1 level and 1 element')
               
        self.vectintegrate = np.vectorize(integrate)
        
    def integrate(self):
        self.vectintegrate(self.sims)
        return

def simapply(config,model,inputdata,**kwargs):
    """Builds an array of models instances using the iterable variables array

    Args:
        config ([dict or path]): base configuration file
        model (cv19model): compartmental model class
        inputdata (cv19data): (optional) Input data for IC and fitting
    """
    def aux(x):
        return(model(config,inputdata,**kwargs,**x))
    aux2 = np.vectorize(aux)
    return(aux2)

def integrate(x):
    """Solves EDOs in models instances

    Args:
        x (cv19model): cv19model instance
    """
    x.integrate()
    return()

def iterate(config, iterables=None,**kwargs,):
    if iterables == None:
        iterables = {}
        nelements = 1
        for key,value in config.items():
            if type(value) == list:
                iterables.update({key:value})
                nelements*=len(value)
        print('Simulating over '+str(len(iterables))+' levels and '+str(nelements)+' elements')        
        
    
    if iterables:
        iterables_aux = deepcopy(iterables)
        iterating = iterables_aux.popitem()
        out = []
        for i in iterating[1]:
            kwargs_aux = deepcopy(kwargs)
            kwargs_aux.update({iterating[0]:i})
            out.append(iterate(config,iterables=iterables_aux,**kwargs_aux))   
        return np.array(out)
    else:
        aux = deepcopy(config)
        for key,value in kwargs.items():
            aux.update({key:value})
        return(aux)

