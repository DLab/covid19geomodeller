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
package_directory = os.path.dirname(os.path.abspath(__file__))

import data.cv19data as cv19data
import utils.cv19functions as cv19functions 
import utils.cv19timeutils as cv19timeutils

"""
# ------------------------------------------------- #   
#                                                   #
#            Configuration File manager             #
#                                                   #
# ------------------------------------------------- #
To Do:
  * Improve this code
  * Improve get_default_parameters. Looks like the current way is deprecated. 
  * Fix meta-population dynamic parameters
   
"""

def loadconfig(sim,config,inputdata,**kwargs):
    # ------------------------------- #
    #        Parameters Load          #
    # ------------------------------- #
    load_default(sim) # Serves for loading optional parameters when they are not defined
    
    # Update parameters using config file
    if type(config) == dict: # If config is a dictionary
        sim.cfg.update(config)
    elif type(config) == str: # If config is a filepath
        sim.cfg.update(toml.load(config))

    # Import static variables
    for key,value in sim.cfg['parameters']['static'].items():
        # Check if the variable is in kwargs
        if key in kwargs:
            value = kwargs[key] # Overwrite by kwargs values
            sim.cfg['parameters']['static'][key]=value # Update config with kwargs values        
        sim.__dict__.update({key:value})
        
    # t_end and t_init are useful when the simulation time should be relative to others
    sim.tsim = sim.t_end - sim.t_init

    # Build dynamic variables
    #TODO: Simplify this code
    if "Metapopulation" in sim.compartmentalmodel:
        # Metapopulation
        for key,value in sim.cfg['parameters']['dynamic'].items():        
            if key in kwargs:
                value = kwargs[key]
                sim.cfg['parameters']['dynamic'][key] = value
            
            if type(value) == list or type(value) == np.ndarray:
                sim.__dict__.update({key:cv19functions.build_metapopulation(value)})
            else:
                sim.__dict__.update({key:cv19functions.build(value)})            
    
    else:
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

    # Convert to datetime format if it's a string
    if sim.initdate and type(sim.initdate) == str:
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
            sim.inputdata.import_data()
            sim.data = sim.inputdata.data
        else: 
            # No data added
            sim.data = None
            sim.inputdata = None

    # ------------------------------- #
    #       Initial conditions        #
    # ------------------------------- #
    sim.initialconditions = sim.cfg['initialconditions']
    for key,value in sim.initialconditions.items():
        # Overwrite by kwargs values
        if key in kwargs:
            value = kwargs[key]
            if type (value) == list:
                value = np.array(value)            
            sim.initialconditions[key] = value
            sim.cfg['initialconditions'][key] = value       

        # Initializing variables with external data if available 
        if type(value) == str:
            # Crear error cuando no haya archivo de datos y fecha inicial
            try:
                sim.__dict__.update({key:sim.data[value][0]})
            except:
                sim.__dict__.update({key:sim.inputdata.__dict__[value]})
        else:
            # Using the values expressed in the configuration file
            if type (value) == list:
                value = np.array(value)
            sim.__dict__.update({key:value})
    return

def load_default(sim):
    sim.cfg = getdefault(sim.compartmentalmodel)

    # Static variables
    for key,value in sim.cfg['parameters']['static'].items():
        sim.__dict__.update({key:value})

    # Dynamic variables
    for key,value in sim.cfg['parameters']['dynamic'].items():          
        sim.__dict__.update({key:cv19functions.build(value)})
        
    
    # Initial conditions        
    sim.initialconditions = sim.cfg['initialconditions']
    for key,value in sim.cfg.items():
        sim.__dict__.update({key:value})    
    return

def getdefault(compartmentalmodel):
    return toml.load(os.path.abspath(package_directory+'/../default_config_files/'+compartmentalmodel+'.toml'))    
    

def saveconfig(filename,config):    
    with open(filename, "w") as toml_file:
        toml.dump(config, toml_file)
    print('Configuration file saved in:\n'+filename)#update
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