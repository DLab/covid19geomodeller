#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import toml
import time
# optimization
import pygmo as pg
#from numpy import linalg as LA

from cv19gm.data.cv19data import ImportData as importdata
from datetime import datetime

from cv19gm.cv19sim import CV19SIM

import cv19gm.utils.cv19functions as cv19functions


def ITWE(sim,data,t=None,rho = 1):
    # Inverse time weighted error
    if type(t) == type(None):
        t = list(range(len(data)))    
    
    err = [((data[i]-sim[t[i]])**2)/(1+t[i])**rho for i in range(len(data))]
    return np.sqrt(np.mean(np.array(err)**2))


def ITWE_log(data,sim,rho = 1):
    # Inverse time weighted error
    err = [((np.log(data[t]+1)-np.log(sim[t]+1))**2)/(1+t)**rho for t in range(len(data))]
    return err

def RMSE(sim,data,t=None,rho = None):
    # rho is just for having the same amount of inputs
    # Inverse time weighted error
    if type(t) == type(None):
        t = list(range(len(data)))    
    
    err = [((data[i]-sim[t[i]])**2)/(1) for i in range(len(data))]
    return np.sqrt(np.mean(np.array(err)**2))

def MSE(sim,data,t=None,rho = None):
    # rho is just for having the same amount of inputs
    # Inverse time weighted error
    if type(t) == type(None):
        t = list(range(len(data)))    
    
    err = [((data[i]-sim[t[i]])**2)/(1) for i in range(len(data))]
    return np.mean(np.array(err)**2)

def P2PE(sim,data,t):
    # Point to point Error
    return [((data[i]-sim[t[i]]))/(1) for i in range(len(data))]


class CV19PARAMFIT():
    def __init__(self, sim, inputdata,verbose = False):
        self.sim = sim


    def ITWErr(sim,data,t=None,rho = 1):
        # Inverse time weighted error
        if type(t) == type(None):
            t = list(range(len(data)))    
        
        err = [((data[i]-sim[t[i]])**2)/(1+t[i])**rho for i in range(len(data))]
        return np.sqrt(np.mean(np.array(err)**2))
    
class BETA_FIT:
    """
    Optimizador de beta con Inverse Time Weighted Error 
    
    
    1. Se deben encontrar las condiciones iniciales para todas las variables
    2. Elegir qué variable se utilizará en el fitness. Voy a partir con I_d
    
    
    """
    
    def __init__(self, I_d, t, cfg, bounds,error = 'RMSE',alpha = 1,**kwargs):
        self.I_d = I_d                  # New daily infected
        self.t = t                    # real time   
        self.bounds = bounds           # Bounds for the variables to optimize
        self.cfg = cfg
        self.alpha=alpha
        
        self.kwargs = kwargs
        
        if error == 'ITWE':
            self.error =  ITWE
        elif error == 'RMSE':
            self.error =  RMSE      
            
        self.ITWE_error = 0
        self.RMSE_error = 0
        self.p2p_error = []
        
        self.mu = mu
        
    def fitness(self,x):
        sim = CV19SIM(self.cfg,beta=x,**self.kwargs)
        sim.solve()
        self.ITWE_error = ITWE(sim.sims[0].I_d,self.I_d,self.t,alpha=self.alpha)
        self.RMSE_error = RMSE(sim.sims[0].I_d,self.I_d,self.t)
        self.p2p_error = P2PE(sim.sims[0].I_d,self.I_d,self.t)    
        return([self.ITWE_error])

    def get_bounds(self):    # mandatory function of Pygmo2
        return(self.bounds)

    def set_bounds(self,bounds):
        self.bounds = bounds
        return(self.bounds)
    
    def gradient(self, x):
        return pg.estimate_gradient_h(lambda x: self.fitness(x), x)


class BETAMU_FIT:
    """
    Optimizador de beta con Inverse Time Weighted Error 
    
    
    1. Se deben encontrar las condiciones iniciales para todas las variables
    2. Elegir qué variable se utilizará en el fitness. Voy a partir con I_d
    
    
    """
    
    def __init__(self, data_I_d, t, cfg, bounds,error='ITWE',rho = 1,**kwargs):
        self.data_I_d = data_I_d                  # New daily infected
        self.t = t                    # real time   
        self.bounds = bounds           # Bounds for the variables to optimize
        self.cfg = cfg
        self.rho=rho
        self.kwargs = kwargs
        
        if error == 'ITWE':
            self.error =  ITWE
        elif error == 'RMSE':
            self.error =  RMSE            
        
    def fitness(self,x):              # mandatory function of Pygmo2
        # We construct the days intervals with the values that are being
        # optimized.

        sim = CV19SIM(self.cfg,beta=x[0],mu=x[1],**self.kwargs)
        sim.solve()
        err = self.error(sim.sims[0].I_d_det,self.data_I_d,t=self.t,rho=self.rho)
        return([err])

    def get_bounds(self):    # mandatory function of Pygmo2
        return(self.bounds)

    def set_bounds(self,bounds):
        self.bounds = bounds
        return(self.bounds)
    
    def gradient(self, x):
        return pg.estimate_gradient_h(lambda x: self.fitness(x), x)