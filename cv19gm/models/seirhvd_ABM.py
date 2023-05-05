import pandas as pd
import numpy as np
from julia import Main as Julia
from julia import Pkg
import cv19gm.utils.cv19files as cv19files

class SEIRHVD_ABM:
    
    def __init__(self, config = None, inputdata = None, verbose = False,  **kwargs):
        self.compartmentalmodel = "SEIRHVD_ABM"
        
        if not config:
            print('Missing configuration file, using default')
    
        if verbose:
            print('Loading configuration file')
                      
        cv19files.loadconfig(self,config,inputdata,**kwargs)
        
        # fill missing atributes
            
        self.nStates = { }
        for state in ['S', 'Sv','E','Ev','I','Im','Ivm','Icr','Ivcr','R','H', 'D']:
            if not hasattr(self, state):
                setattr(self, state, 0)
            self.nStates[state] = int(getattr(self, state))
            
        if not hasattr(self, 'beta_nn') or not hasattr(self, 'beta_nv') or not hasattr(self, 'beta_vn') or not hasattr(self, 'beta_vv'):
            self.beta_nn = self.beta
            self.beta_nv = self.beta
            self.beta_vn = self.beta
            self.beta_vv = self.beta
            
        self.chanceInfect = {
            (False, False): self.beta_nn,
            (False, True): self.beta_nv,
            (True, False): self.beta_vn,
            (True, True): self.beta_vv
        }
        
        self.tRecover = {
            True: self.tImv_R(0),
            False: self.tIm_R(0)
        }
        self
        
        self.vRecover = {
            True: self.vImv_R(0),
            False: self.vIm_R(0)
        }
        
        self.tCriticalDie = {
            True: self.tIv_D(0),
            False: self.tIcr_D(0)
        }
        
        self.vCriticalDie = {
            True: self.vIv_D(0),
            False: self.vIcr_D(0)
        }
             
    def run(self):
        Pkg.activate("./julia")
        Julia.include("./julia/run.jl")
        
        data = Julia.run_SEIRHVD(
            self.nStates,
            
            self.alpha,
            self.chanceInfect,
            self.pIcr_H,
            self.pH_D,
            self.vac_d,
            
            self.tE_I(0),
            self.vE_I(0),
            self.tRecover,
            self.vRecover,
            self.tCriticalDie,
            self.vCriticalDie,
            self.tH_D(0),
            self.vH_D(0),
            self.tH_R(0),
            self.vH_R(0),
            self.tR_S(0),
            self.vR_S(0),
            
            self.stepsPerDay,
            self.t_end - self.t_init,
            
            isGraphSpace = self.network,
            startDay = self.t_init + 1)
        
        #set attributes
        for state in ['S', 'E', 'Im', 'Icr', 'R', 'H', 'D']:
            setattr(self, state, dict())
            getattr(self, state)[False] = data['totals'][state, False]
            getattr(self, state)[True] = data['totals'][state, True]
            
            dailyname = state + '_d'
            setattr(self, dailyname, dict())
            getattr(self, dailyname)[False] = data['daily'][state, False]
            getattr(self, dailyname)[True] = data['daily'][state, True]
        
        
        #set totals
        self.totals = pd.DataFrame()
        for (sym, vac) in data['totals'].keys():
            name = sym + ('v' if vac else '')
            self.totals[name] = data['totals'][sym,vac]
            #dic = np.array(data['totals'][sym,vac])
            #setattr(self, name, dic)
            
        #set daily
        self.daily = pd.DataFrame()
        for (sym,vac) in data['daily'].keys():
            name = sym + ('v' if vac else '') + '_d'
            self.daily[name + '_d'] = data['daily'][sym,vac]
            #dic = dict()
            #dic[sym,vac] = np.array(data['daily'][sym,vac])
            #setattr(self, name, dic)
            
        self.results = pd.concat([self.totals, self.daily], axis=1)
        
        

        
        
        
        
        
    
    
        

        
        
        
        
        
    