import numpy as np
import pandas as pd
import toml

import cv19gm.utils.cv19files as cv19files
from copy import deepcopy



"""
CV19SIM interphaces the user with all the cv19gm library tools.

Todo: [ ] Build resume function
        * tipo de modelo
        * variables a iterar
        * Condiciones iniciales reales 
Todo: [ ] Build a new module for multiple simulations. 
Todo: [ ] Add metapopulation models
Todo: [ ] Make variables accessible from parent object
Todo: [ ] simplificar la vectorización de la función de integración

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
            
        elif compartmentalmodel == 'SEIR_META':    
            self.modelname = compartmentalmodel
            from cv19gm.models.seir_meta import SEIR_META
            compartmentalmodel = SEIR_META
            
        else:
            raise('Incorrect model')

        # I need to overload the current class with the following instance of the model       
        self.sims = compartmentalmodel(config,inputdata,**kwargs)
        
      
        
    def solve(self):
        self.sims.solve()
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
        
