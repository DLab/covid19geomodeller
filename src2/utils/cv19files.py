#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.special import expit
import json
import pandas as pd
import toml

import os
import sys
path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(1, path)

import data.cv19data as cv19data
import utils.cv19functions as cv19functions 
import utils.cv19timeutils as cv19timeutils


"""
# ------------------------------------------------- #   
#                                                   #
#                   File manager                    #
#                                                   #
# ------------------------------------------------- #
Todo:
    [] Capacidad de guardar un archivo de configuración de un objeto de simulación en un archivo toml
Functions:
    loadconfig: Load parameters from a configuration file or object
    saveconfig: Save parameters into a configuration file
    loaddata: Load data from file
    savedata: Save data to file

 
"""
def loadconfig(sim,config,inputdata,**kwargs):
    """[summary]

    Args:
        sim ([type]): [description]
        config ([type]): [description]
        inputdata ([type]): [description]
    """
    # ------------------------------- #
    #        Parameters Load          #
    # ------------------------------- #
    if type(config) == dict:
        sim.cfg = config    
    else:
        sim.cfg = toml.load(config)
            
    # sim
    sim.model = sim.cfg['model']

    # Import fixed variables
    for key,value in sim.cfg['parameters']['static'].items():
        if key in kwargs:
            value = kwargs[key]
            sim.cfg['parameters']['static'][key]=value
        sim.__dict__.update({key:value})
        
    
    sim.tsim = sim.t_end - sim.t_init

    # Build functions
    for key,value in sim.cfg['parameters']['dynamic'].items():        
        if key in kwargs:
            value = kwargs[key]
            sim.cfg['parameters']['dynamic'][key] = value                
        sim.__dict__.update({key:cv19functions.build(value)})
        
    # Ephemeris
    if 'ephemeris' in sim.cfg:
        sim.ephemeris = sim.cfg['ephemeris']

    # Data:
    for key,value in sim.cfg['data'].items():
        if key in kwargs:
            value = kwargs[key]
            sim.cfg['data'][key] = value        
        sim.__dict__.update({key:value})

    if sim.initdate:
        sim.initdate = cv19timeutils.txt2Datetime(sim.initdate) 

    # cv19data object
    if inputdata:
        sim.inputdata = inputdata
        sim.data = sim.inputdata.data
        sim.initdate = sim.inputdata.initdate
        #if not type(sim.initdate) == datetime.datetime:
        #    sim.initdate = cv19timeutils.txt2Datetime(sim.initdate)                

        #sim.country = sim.inputdata.country 
        sim.state = sim.inputdata.tstate 
        #sim.county = sim.inputdata.county             
        #if sim.inputdata.loc_name:
            #sim.loc_name = sim.inputdata.loc_name

    
    else:
        # Local file
        if sim.datafile:
            # Falta construir funcion para cargar datos desde archivos, pero primero tengo que construir ese archivo 
            sim.data = pd.DataFrame(sim.cfg['data']['datafile'])
            sim.inputdata = None
        elif sim.importdata:                	            
            #sim.inputdata = cv19data.ImportData(country=sim.country,state=sim.state,county=sim.county,healthservice=sim.healthservice,initdate=sim.initdate,user = None,password = None)
            sim.inputdata = cv19data.ImportData(tstate=sim.state,initdate=sim.initdate,user = None,password = None)
            sim.inputdata.importdata()
            sim.data = sim.inputdata.data
        else: 
            #print('No external data added')
            sim.data = None
            sim.inputdata = None

    # ------------------------------- #
    #       Initial conditions        #
    # ------------------------------- #
    sim.initialconditions = sim.cfg['initialconditions']
    for key,value in sim.initialconditions.items():
        if key in kwargs:
            value = kwargs[key]
            sim.initialconditions[key] = value
            sim.cfg['initialconditions'][key] = value
        if type(value) == str:
            # Crear error cuando no haya archivo de datos y fecha inicial
            sim.__dict__.update({key:sim.data[value][0]})
        else:
            sim.__dict__.update({key:value})
    return

def saveconfig(filename,config):    
    with open(filename, "w") as toml_file:
        toml.dump(config, toml_file)
    print('Configuration file saved in:\n'+filename)#update
    return 

def loaddata(filename):
    return

def savedata(data=None,simulation=None,filename=None):
    if data:
        pass
    elif simulation:
        pass
    else:
        print("No data given")
        # raise error
        return
    return 


def unwrapconfig(config,**kwargs):
    # ------------------------------- #
    #        Parameters Load          #
    # ------------------------------- #
    out = {}

    if type(config) == dict:
        cfg = config    
    else:
        cfg = toml.load(config)
            
    # Model
    out.update({'model':cfg['model']})

    # Import fixed variables
    for key,value in cfg['parameters']['static'].items():
        if key in kwargs:
            value = kwargs[key]
            cfg['parameters']['static'][key]=value
        out.update({key:value})
        
    

    # Build functions
    for key,value in cfg['parameters']['dynamic'].items():        
        if key in kwargs:
            value = kwargs[key]
            cfg['parameters']['dynamic'][key] = value                
        out.update({key:value})
        
    # Ephemeris
    if 'ephemeris' in cfg:
        out.update({'ephemeris':cfg['ephemeris']})

    # Data:
    for key,value in cfg['data'].items():
        if key in kwargs:
            value = kwargs[key]
            cfg['data'][key] = value        
        out.update({key:value})

    

    # ------------------------------- #
    #       Initial conditions        #
    # ------------------------------- #    
    for key,value in cfg['initialconditions'].items():
        if key in kwargs:
            value = kwargs[key]            
            cfg['initialconditions'][key] = value
        out.update({key:value})
    return out