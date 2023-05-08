import numpy as np
import pandas as pd
import toml

import cv19gm.utils.cv19files as cv19files
from copy import deepcopy



"""
CV19SIM interphaces the user with all the cv19gm library tools.
It also adds the capability of performing multiple simulations in order to study the behavior of some parameters.  

Todo: [ ] Build resume function
        * tipo de modelo
        * variables a iterar
        * Condiciones iniciales reales 
Todo: [ ] Build a new module for multiple simulations. 
Todo: [ ] Make variables accessible from parent object
Todo: [ ] simplificar la vectorización de la función de integración
Todo: [ ] Paralelizar la integración de las EDOs dentro de lo posible
Todo: [ ] Solve está siendo aplicado más de una vez. Ver donde ocurre eso! 

"""

class CV19SIM():
    def __init__(self,config = None, inputdata = None, compartmentalmodel=None, verbose=False,**kwargs):
        
        # Leer el archivo de configuracion
        if not config:
            if not compartmentalmodel:
                raise Exception("Missing compartmental model definition")
            print("Using default configuration file for "+compartmentalmodel+" model")
            config = cv19files.getdefault(compartmentalmodel)
        else:
            if not type(config) == dict:
                config = toml.load(config)
            else:
                config = deepcopy(config)
            
        aux = cv19files.unwrapconfig(config,**kwargs)        
        if not compartmentalmodel:
            compartmentalmodel = aux['model']['model']
        
        if compartmentalmodel == 'SEIR':
            self.modelname = compartmentalmodel
            from cv19gm.models.seir import SEIR
            compartmentalmodel = SEIR
             
        elif compartmentalmodel == 'SEIRHVD':
            self.modelname = compartmentalmodel
            from cv19gm.models.seirhvd import SEIRHVD
            compartmentalmodel = SEIRHVD

        elif compartmentalmodel == 'SIR':
            self.modelname = compartmentalmodel
            from cv19gm.models.sir import SIR
            compartmentalmodel = SIR

        elif compartmentalmodel == 'SEIRTQ':
            self.modelname = compartmentalmodel
            from cv19gm.models.seirtq import SEIRTQ
            compartmentalmodel = SEIRTQ
        else:
            raise('Incorrect model')


        sims = []
        self.iterables = {}
        
        # Find iterable parameters:
        for key,value in aux.items():
            if type(value)==list:
                if verbose:
                    print(key+':'+str(value))
                self.iterables.update({key:value})

        #print('There are '+ str(len(iterables))+' iterable parameters')

        if self.iterables:
            # Pop iterables from kwargs
            for key in self.iterables:
                if key in kwargs:
                    kwargs.pop(key)
            expandediterables = iterate(self.iterables,verbose=verbose)
            buildsim = simapply(config,compartmentalmodel,inputdata,**kwargs)            
            self.sims = buildsim(expandediterables)
        else:
            self.sims = [compartmentalmodel(config,inputdata,**kwargs)]
            if verbose:
                print('Simulating over 1 level and 1 element')
               
        self.vectsolve = np.vectorize(solve)
        
        if verbose:
            print(str(np.prod(np.shape(self.sims)))+" models created")
        
    def integrate(self):
        print('The use of integrate() is now deprecated. Use solve() instead.')
        self.vectsolve(self.sims)
        return        
        
    def solve(self):
        self.vectsolve(self.sims)
        return

    def resume(self):
        """Resume:
        Prints a resume of the object with

        Model type:
        Simulated:
        Number of simulations: 
        Iterated variables:
        RealData:
            * CUT
            * Dates
        """
        print('Model: '+self.modelname)
        print('Iterables: '+self.iterables)
        return
    
    def extract(self):
        """Extract variables from simulation objects
        """
        shape = np.shape(self.sims)
        
    def expose_variable(self,variable):
        """Expose variables so they can be accessed directly from the main class object

        Args:
            variable (string): Variable to be exposed
        """
        setattr(self,variable, np.reshape(list(map(lambda sim: sim.__dict__[variable],self.sims.flatten())),np.shape(self.sims)))
        return
        
def simapply(config,compartmentalmodel,inputdata,**kwargs):
    """Builds an array of models instances using the iterable variables array

    Args:
        config ([dict or path]): base configuration file
        model (cv19model): compartmental model class
        inputdata (cv19data): (optional) Input data for IC and fitting
    """
    def aux(x):
        return(compartmentalmodel(config,inputdata,**kwargs,**x))
    aux2 = np.vectorize(aux)
    return(aux2)

def solve(x):
    """Solves EDOs in models instances

    Args:
        x (cv19model): cv19model instance
    """
    x.solve()
    return()

def iterate(config, iterables=None,verbose=False,**kwargs,):
    if iterables == None:
        iterables = {}
        nelements = 1
        for key,value in config.items():
            if type(value) == list:
                iterables.update({key:value})
                nelements*=len(value)
        if verbose:
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


