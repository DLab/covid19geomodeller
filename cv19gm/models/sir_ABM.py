import pandas as pd
from julia import Main as Julia
from julia import Pkg

import cv19gm.utils.cv19files as cv19files

class SIR_ABM:
    
    def __init__(self, config = None, inputdata = None, verbose = False,  **kwargs):
        self.compartmentalmodel = "SIR_ABM"
        
        if not config:
            print('Missing configuration file, using default')

        if verbose:
            print('Loading configuration file')
                      
        cv19files.loadconfig(self,config,inputdata,**kwargs)
        
        
        # fill missing atributes

        if not hasattr(self, 'S') or not self.S:
            self.S = self.population - self.I - self.R
            
        self.S = int(self.S)
        self.I = int(self.I)
        self.R = int(self.R)
        
        self.days = self.t_end - self.t_init if hasattr(self, 't_end') else None
        
    def run(self):
        Julia.using("Pkg")
        Pkg.activate("./julia")
        Julia.include("./julia/run.jl")
        
        data = Julia.run_SIR(self.S, self.I, self.R, 
                             self.alpha, self.beta, self.tI_R(0), self.vI_R(0), 5, 
                             isGraphSpace=self.network,
                             days = self.days, startDay = self.t_init + 1)
        
        self.daily = data['daily'].rename(columns = {'S': 'S_d', 'I': 'I_d', 'R': 'R_d'})
        self.totals = data['totals']
        self.results = pd.concat([self.totals,self.daily], axis=1)
        
        for col in self.results:
            setattr(self, col, self.results[col].values)