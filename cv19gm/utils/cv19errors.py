#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
"""
# -------------------------------- #
#     Error measuring functions    #
# -------------------------------- #
"""
def cv19errorbuild(input):
    """
    Function builder
    # crear un iterador que recorra el input y cree la función a partir de un comando exec:
    # Acepta diccionarios o strings con forma de diccionario
    """
    if callable(input):        
        return input
    
    elif type(input)==str:        
        print('input dict')
        return locals()[input]


# -------------------------------- #
#          Global Errors           #
# -------------------------------- #

def ITWE(sim, data, t_data=None, rho=1):
    # Inverse time weighted error
    if type(t_data) == type(None):
        t_data = list(range(len(data)))
    err = [((data[i] - sim[t_data[i]]) ** 2) / (1 + t_data[i]) ** rho for i in range(len(data))]
    return np.sqrt(np.mean(np.array(err) ** 2))


def ITWE_log(sim, data, t_data=None, rho=1):
    # Log Inverse time weighted error
    if type(t_data) == type(None):
        t_data = list(range(len(data)))
    err = [((np.log(data[i]+1) - np.log(sim[t_data[i]]+1)) ** 2) / (1 + t_data[i]) ** rho for i in range(len(data))]
    return np.sqrt(np.mean(np.array(err) ** 2))

def RMSE(sim, data, t_data=None, rho=None):
    # Root Mean Squared Error
    if type(t_data) == type(None):
        t_data = list(range(len(data)))

    err = [((data[i] - sim[t_data[i]]) ** 2) / (1) for i in range(len(data))]
    return np.sqrt(np.mean(np.array(err) ** 2))

def MSE(sim, data, t_data=None, ):
    # Mean Square Error
    if type(t_data) == type(None):
        t_data = list(range(len(data)))

    err = [((data[i] - sim[t_data[i]]) ** 2) / (1) for i in range(len(data))]
    return np.mean(np.array(err) ** 2)

def RRMSE(sim, data, t_data=None):
    # Relative Root Mean Square Error
    if type(t_data) == type(None):
        t_data = list(range(len(data)))

    err = [((data[i] - sim[t_data[i]]) ** 2) / (1 + data[i]) for i in range(len(data))]
    return np.sqrt(np.mean(np.array(err) ** 2))


# -------------------------------- #
#        Residual Metrics          #
# -------------------------------- #

def LAE(sim, data, t_data=None):
    # Local Absolute Error. Output is a vector
    if type(t_data) == type(None):
        t_data = list(range(len(data)))    
    return [np.abs(data[i] - sim[t_data[i]])  for i in range(len(data))]

def LRAE(sim, data, t_data=None):
    # Local Relative Error. Output is a vector
    if type(t_data) == type(None):
        t_data = list(range(len(data)))
    return [(np.abs(data[i] - sim[t_data[i]])) / (1 + data[i]) for i in range(len(data))]    

def LCE(sim, data, t_data=None):
    """Local Cummulative Error 
    To do
    Args:
        sim (list): Simulation data
        data (list): Data to compare
        t_data (list, optional): _description_. Defaults to None.

    Returns:
        list: accumulative error
    """
    # Local Integrative Error. Output is a vector
    if type(t_data) == type(None):
        t_data = list(range(len(data)))
    return np.cumsum([data[i] - sim[t_data[i]] for i in range(len(data))])
    
def LRCE(sim, data, t_data=None):
    """Local Relative Cummulative Error 
    To do
    Args:
        sim (list): Simulation data
        data (list): Data to compare
        t_data (list, optional): _description_. Defaults to None.

    Returns:
        list: accumulative error
    """
    # Local Integrative Error. Output is a vector
    if type(t_data) == type(None):
        t_data = list(range(len(data)))
    return np.cumsum([(data[i] - sim[t_data[i]]) / (1 + data[i]) for i in range(len(data))])
    