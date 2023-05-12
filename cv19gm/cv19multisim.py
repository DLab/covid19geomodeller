import numpy as np
import toml

import cv19gm.utils.cv19files as cv19files
from copy import deepcopy

"""
CV19MULTISIM 
Performs multiple simulations for sensitivity analysis on single population models.


Todo: [ ] Simplify the code
Todo: [ ] Build resume function
Todo: [ ] Make variables accessible from parent object
Todo: [ ] Add metapopulation models (?)
"""


class CV19MULTISIM():
    def __init__(self,config = None, compartmentalmodel=None, verbose=False,**kwargs):        
        # Read the configuration file
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
            compartmentalmodel = aux['model']['name']
        
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
            raise('Incorrect model: '+str(compartmentalmodel))

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
            buildsim = simapply(config,compartmentalmodel,**kwargs)            
            self.sims = buildsim(expandediterables)
        else:
            self.sims = [compartmentalmodel(config,**kwargs)]
            if verbose:
                print('Simulating over 1 level and 1 element')
        
        if verbose:
            print(str(np.prod(np.shape(self.sims)))+" models created")
        
    def solve(self):
        self.vectsolve = np.vectorize(solve)
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
        
    def expose_variable(self,variable):
        """Expose variables so they can be accessed directly from the main class object

        Args:
            variable (string): Variable to be exposed
        """
        setattr(self,variable, np.reshape(list(map(lambda sim: sim.__dict__[variable],self.sims.flatten())),np.shape(self.sims)))
        return
        
def simapply(config,compartmentalmodel,**kwargs):
    """Builds an array of models instances using the iterable variables array

    Args:
        config ([dict or path]): base configuration file
        model (cv19model): compartmental model class
    """
    def aux(x):
        return(compartmentalmodel(config,**kwargs,**x))
    #aux2 = np.vectorize(aux)
    return(np.vectorize(aux))

def solve(x):
    """Solves EDOs in models instances

    Args:
        x (cv19model): cv19model instance
    """
    x.solve()
    return()


def iterate(config, iterables=None,verbose=False,**kwargs,):
    """Recursive function for iterating over the objective parameters

    Args:
        config (cv19config): cv19config instance
        iterables (list, optional): list of parameters being iterated over. Defaults to None.
        verbose (bool, optional): Verbose. Defaults to False.

    Returns:
        np.array: Return array of cv19sim instances
    """
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
