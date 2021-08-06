#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.special import expit
import json


"""
# ------------------------------------------------- #   
#                                                   #
#        Time dependent function constructors       #
#                                                   #
# ------------------------------------------------- #

To Do:    
    * Agregar tipos de subidas y bajadas al Events function
    * Reducir las funciones de subida y bajada continuas en 1 con tipo de subida y bajada
    * Revisar el valor por defecto de events
    * Revisar si se puede simplificar la llamada a events 
    * Considerar el crear funciones directamente, sin pasar por ese formato que cree   
"""
def build(input):
    """
    Function builder
    # crear un iterador que recorra el input y cree la función a partir de un comando exec:
    # Acepta diccionarios o strings con forma de diccionario
    """
    
    if type(input)==str:
        input_dict = json.loads(input)
    elif type(input)==dict:
        input_dict = input.copy()
    elif callable(input):
        setattr(locals()['out'],'constructor','User defined Function')
        return input
    elif type(input) == tuple:
        return addbuild(*input)
    elif type(input) == list:
        # For building iterable simulations
        return input        
    else:        
        def out(t):
            return input
        setattr(locals()['out'],'constructor',str(input))            
        return out
    
    #try:
    #    print("Executing "+input_dict['function'])        
    #except:
    #    raise SyntaxError("No function defined")
    aux = 'out=' + input_dict['function']+"("
    del input_dict['function']
    for key, value in input_dict.items():
        aux +=key+"="+str(value)+',' 
    aux +=')'
    #print(aux)
    ldict={}
    exec(aux,globals(),ldict)
    out = ldict['out']
    setattr(locals()['out'],'constructor',str(input))
    return locals()['out']

# Build + function addition.
def addbuild(*input):
    """
    Creates a function Build a functions and adds t 
    # Construir funcion que sume multiples funciones 
    """
    aux = []
    for i in input:
        aux.append(build(i))
    return functionAddition(*aux)

# Function addition
def functionAddition(*functions):
    """
    Creates the function that results adding the functions received. 
    Receives a list where each element is a function that receives 1 argument
    """
    def f(t):
        aux = 0
        for i in functions:    
            aux += i(t)
        return aux
    return f


def polyfit(values,time = None,degree=4,endvalue_index = -5):
    """
    polyfit:
    Function data fits real data with a polynom of a given degree, and then projects the future values with it.
    values: values to fit
    time: time array from data to fit. If no time array given, daily data is assumed 
    degree: polynom degree    
    """
    if not time:
        time = np.array(range(len(values)))
    datamodel = np.poly1d(np.polyfit(time, values, degree))
    # transition time from data to fixed value
    tchange = time[-1]#[tchange]
    endvalue = np.mean(values[endvalue_index:])
    f_out=lambda t: datamodel(t)*(1-expit(t-tchange)) + expit(t-tchange)*endvalue
    return f_out

# TODO: Add default value
# TODO: Change functions to *args

def Events(values,days,functions = []):
    """
    Event creator function. Create a time dependent function that returns the values setted in the 
    values vector for the periods specified in the days vector. 


    Args:
        values (list): list with the values for the different intervals
        days (list): ist of lists of len == 2 with the values for the function on that v.
        functions (list, optional): [description]. Defaults to [].

    Returns:
        [type]: [description]
    """

    # for one interval
    #if type(not days[0]) == list and len(days == 2):
    #    days = [[days[0],days[1]]]

    # define default function
    f = lambda t:0

    aux_f = [f]
    for i in functions:
        aux_f.append(i)
        

    
    for i in range(len(values)):
        def auxf(t,j=i):
            return values[j]*(expit(10*(t-days[j][0])) - expit(10*(t-days[j][1]))) 
        aux_f.append(auxf)

    out = functionAddition(*aux_f)    
    return out


    """
    days = [[3,10],[15,25],[8,17]]
    values = [1,5,3]
    t = np.array(np.arange(0,30,0.1))
    plt.plot(t,sumfunct(t))
    plt.show()
    """

# Periodic square function
def square(min_val=0,max_val=1,period=14,init=0,end=1000,off_val=0,initphase='min',duty=0.5):
    # phase en fraccion del periodo
    def f(t): 
        return signal.square(t,duty)
    
    if initphase == 'min':
        phi = np.pi*(2*init/period - 1.1)
    elif initphase == 'max':
        phi = 2*np.pi*init/period
    else:
        phi = 2*np.pi*(init/period - initphase)
    def aux(t):    
        return (expit(10*(t-init)) - expit(10*(t-end)))*((max_val-min_val)/2*(f(2*np.pi / period * t - phi))+(max_val+min_val)/2) + \
        (1-(expit(10*(t-init)) - expit(10*(t-end))))*off_val

    return aux

# Sine function
def sine(min_val=0,max_val=1,period=14,init=0,end=1000,off_val=0,initphase='min'):
    # phase en fraccion del periodo
    def f(t): 
        return np.cos(t)
    
    if initphase == 'min':
        phi = np.pi*(2*init/period - 1)
    elif initphase == 'max':
        phi = 2*np.pi*init/period
    else:
        phi = 2*np.pi*(init/period - initphase)
    def aux(t):    
        return (expit(10*(t-init)) - expit(10*(t-end)))*((max_val-min_val)/2*(f(2*np.pi / period * t - phi))+(max_val+min_val)/2) + \
        (1-(expit(10*(t-init)) - expit(10*(t-end))))*off_val

    return aux

# Sawtooth function
def sawtooth(min_val=0,max_val=1,period=14,init=0,end=1000,off_val=0,initphase='min',width=1):
    # phase en fraccion del periodo
    def f(t): 
        return signal.sawtooth(t,width)
    
    if initphase == 'min':
        phi = np.pi*(2*init/period - 0.5/period)
    elif initphase == 'max':
        phi = 2*np.pi*init/period
    else:
        phi = 2*np.pi*(init/period - initphase)
    def aux(t):    
        return (expit(10*(t-init)) - expit(10*(t-end)))*((max_val-min_val)/2*(f(2*np.pi / period * t - phi))+(max_val+min_val)/2) + \
        (1-(expit(10*(t-init)) - expit(10*(t-end))))*off_val

    return aux

# Increasing functions
def increase_linear(t0,t1,t2,maxvalue = 0,increaserate=1):
    """
    Creates a function that starts growing linearly from t = t0, with the specified increase rate until it reaches the final
    value, which is either maxvalue or increaserate*(t1-t0). After this it keeps this value until t2, and then goes back to 0
    """
    if maxvalue == 0:
        #a = np.polyfit([t0,t0+1],[0,increaserate],1)
        maxvalue = increaserate*(t1-t0)
    
    a = np.polyfit([t0,t1],[0,maxvalue],1)        
    f = lambda t: np.poly1d(a)(t)

    aux = lambda t: f(t)*(expit(10*(t-t0)) - expit(10*(t-t1))) + maxvalue*(expit(10*(t-t1)) - expit(10*(t-t2)))
    return aux

def increase_quadratic(t0,t1,t2,maxvalue = 1):
    """
    Creates a function that starts growing quadratically from t = t0, until it reaches the maxvalue. 
    After this it keeps this value until t2, and then goes back to 0
    """
    a = np.polyfit([t0,t1,2*t0-t1],[0,maxvalue,maxvalue],2)
    f = lambda t: np.poly1d(a)(t)

    aux = lambda t: f(t)*(expit(10*(t-t0)) - expit(10*(t-t1))) + maxvalue*(expit(10*(t-t1)) - expit(10*(t-t2)))
    return aux

def increase_sigmoid(t0,t1,t2,maxvalue = 1):
    """
    Creates a function that grows like a sigmoid from t = t0, until it reaches the maxvalue. 
    After this it keeps this value until t2, and then goes back to 0
    """                  
    def f(t):
        return maxvalue*(expit((t-t0-4)*8/(t1-t0)) - expit(8*(t-t2)))      
    return f

# Saturation function
def saturation(satfunct,gw = 20):
    """
    Saturation function
    Binary function that indicates when the sum of functions are bigger than the saturation function.
    input: 
        satfunct: Upper limit function
    return:
        0 when the functions given are smaller than the saturation function 
        1 when they are bigger or equal
    """
    
    if not callable(satfunct):    
        satfunct = build(satfunct)
    def aux(t,*functions):
        f = functionAddition(*functions)
        return(expit(gw*(f-satfunct(t))))
    return aux
    