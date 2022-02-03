#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.special import expit
import json
import pandas as pd
import ast

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
        #input_dict = json.loads(input)
        print('input dict')
        input_dict = ast.literal_eval(input)
    elif type(input)==dict:
        input_dict = input.copy()
    elif callable(input):
        out = input
        setattr(locals()['out'],'constructor','User defined Function')
        return input
    elif type(input) == tuple:
        return build_add(*input)
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
def build_add(*input):
    """
    Creates a function Build a functions and adds t 
    # Construir funcion que sume multiples funciones 
    """
    aux = []
    for i in input:
        aux.append(build(i))
    return func_add(*aux)

# Function addition
def func_add(*functions):
    """
    Creates the function that results adding the functions received. 
    Receives a list where each element is a function that receives 1 argument
    """
    initfunctions = []
    for i in functions:
        initfunctions.append(build(i)) 

    def f(t):
        aux = 0
        for i in initfunctions:
            aux += i(t)
        return aux
    return f


def events(values,days,default=0,*functions):
    """Event creator function. Create a time dependent function that returns the values (or functions) given in the 
    values vector for the intervals specified in the days vector. 


    Args:
        values (list): list with the values for the different intervals
        days (list): list of lists of len == 2 with the values for the function on that v.
        default (float, optional): value that the function takes in undefined intervals. Defaults to 0
        
        functions: extra functions to add with no interval definition

    Returns:
        out (function): function created
    """

    gw = 20

    f = lambda t:0

    aux_f = [f]
    for i in functions:
        aux_f.append(i)
        

    fvalues = []
    for i in values:
        fvalues.append(build(i))

    for i in range(len(values)):
        def auxf(t,j=i):
            return (fvalues[j])(t)*(expit(gw*(t-days[j][0])) - expit(gw*(t-days[j][1]))) 
        aux_f.append(auxf)

    #
    # Default value
    #

    # find overlapping intervals
    dayssort = np.array(days)    
    dayssort = dayssort[np.argsort(dayssort[:,0])]

    overlapdays = [dayssort[0,0]]
    aux = dayssort[0,1]
    for i in range(1,len(dayssort)):
        for j in range(2):
            if not j and aux < dayssort[i,j]:
                overlapdays.append(aux)
                overlapdays.append(dayssort[i,j])
                aux = dayssort[i,1]
                break
            elif j and aux < dayssort[i,j]:
                aux = dayssort[i,j]
    overlapdays.append(aux)
    overlapdays = np.reshape(overlapdays,(int(len(overlapdays)/2),2))
    
    # Create default function
    for i in range(len(overlapdays)):
        def auxf(t,j=i):
            return -default*((expit(gw*(t-overlapdays[j][0])) - expit(gw*(t-overlapdays[j][1])))) 
        aux_f.append(auxf)
    aux_f.append(lambda t:default)
    
    # Create final function adding them all
    out = func_add(*aux_f)    
    return out



# Periodic square function
def square(min_val=0,max_val=1,period=14,t_init=0,t_end=1000,default=0,initphase=0,duty=0.5):
    # phase en fraccion del periodo
    def f(t): 
        return signal.square(t,duty)
    

    # Start with max value
    if initphase:
        phi = 2*np.pi*t_init/period
    # Start with min value
    else:
        phi = np.pi*(2*t_init/period - 1.1)    
    
    
    def aux(t):    
        return (expit(10*(t-t_init)) - expit(10*(t-t_end)))*((max_val-min_val)/2*(f(2*np.pi / period * t - phi))+(max_val+min_val)/2) + \
        (1-(expit(10*(t-t_init)) - expit(10*(t-t_end))))*default

    return aux

# Sine function
def sine(min_val=0,max_val=1,period=14,t_init=0,t_end=1000,default=0,initphase=0):
    # phase en fraccion del periodo
    def f(t): 
        return np.cos(t)
    
    # Start with max value
    if initphase:
        phi = 2*np.pi*t_init/period 
    # Start with min value   
    else:
        phi = np.pi*(2*t_init/period - 1)

    def aux(t):    
        return (expit(10*(t-t_init)) - expit(10*(t-t_end)))*((max_val-min_val)/2*(f(2*np.pi / period * t - phi))+(max_val+min_val)/2) + \
        (1-(expit(10*(t-t_init)) - expit(10*(t-t_end))))*default

    return aux

# Sawtooth function
def sawtooth(min_val=0,max_val=1,period=14,t_init=0,t_end=1000,default=0,initphase=0,width=1):
    # phase en fraccion del periodo
    def f(t): 
        return signal.sawtooth(t,width)
    
    # Start with max value
    if initphase:
        phi = 2*np.pi*t_init/period
    # Start with max value
    else:
        phi = np.pi*(2*t_init/period - 0.5/period)

    def aux(t):    
        return (expit(10*(t-t_init)) - expit(10*(t-t_end)))*((max_val-min_val)/2*(f(2*np.pi / period * t - phi))+(max_val+min_val)/2) + \
        (1-(expit(10*(t-t_init)) - expit(10*(t-t_end))))*default

    return aux

# Value transition functions
def linear_transition(t_init,t_end,initvalue=0,endvalue = 1):
    """linearTransition
    Creates a function which performs a linear transition from initvalue to endvalue between t_init and t_end.
    Args:
        t_init (int): Transition beginning
        t_end (int): Transition end
        initvalue (int, optional): Initial value. Defaults to 0.
        endvalue (int, optional): End value. Defaults to 1.

    Returns:
        function: function that performs the linear transition
    """
    a = np.polyfit([t_init,t_end],[initvalue,endvalue],1)        
    f = lambda t: np.poly1d(a)(t)
    out = lambda t: initvalue*expit(10*(t_init-t)) + f(t)*(expit(10*(t-t_init)) - expit(10*(t-t_end))) + endvalue*expit(10*(t-t_end)) 
    return out


def quadratic_transition(t_init,t_end,initvalue=0,endvalue = 1,concavity=0):
    """quadraticTransition
    Creates a function which performs a quadratic transition from initvalue to endvalue between t_init and t_end.

    Args:
        t_init (int): Transition beginning
        t_end (int): Transition end
        initvalue (int, optional): Initial value. Defaults to 0.
        endvalue (int, optional): End value. Defaults to 1.
        concavity (int, optional): Function's concavity. 0: convex, 1: concave. Defaults to convex.

    Returns:
        function: function that performs the quadratic transition
    """

    # convex increase or concave decrease
    if (not concavity and endvalue>=initvalue) or (concavity and initvalue>endvalue):
        t_aux = 2*t_init - t_end 
        v_aux = endvalue
        
    # concave increase or convex decrease
    elif (concavity and endvalue>=initvalue) or (not concavity and initvalue>endvalue):
        t_aux = 2*t_end - t_init 
        v_aux = initvalue
        
    a = np.polyfit([t_init,t_end,t_aux],[initvalue,endvalue,v_aux],2)
    f = lambda t: np.poly1d(a)(t)
    
    out = lambda t: initvalue*expit(10*(t_init-t)) + f(t)*(expit(10*(t-t_init)) - expit(10*(t-t_end))) + endvalue*expit(10*(t-t_end)) 
    return out    


def sigmoidal_transition(t_init,t_end,initvalue=0,endvalue = 1,gw=8):
    """Sigmoidal Transition
    Creates a function which performs a sigmoidal transition from initvalue to endvalue between t_init and t_end.

    Args:
        t_init (int): Transition beginning
        t_end (int): Transition end
        initvalue (int, optional): Initial value. Defaults to 0.
        endvalue (int, optional): End value. Defaults to 1.
        gw (float, optional): Gain weight. Calibrates the "strength" of the sigmoid change 

    Returns:
        function: function that performs the sigmoidal transition
    """

    out = lambda t:  initvalue + (endvalue-initvalue)*(expit((t-(t_init+t_end)/2)*gw/(t_end-t_init))) 
    return out 


# Value transition functions
def transition(t_init,t_end,type = 'linear', initvalue=0,endvalue = 1, concavity=0, gw=8):
    """linearTransition
    Creates a function which performs a linear transition from initvalue to endvalue between t_init and t_end.
    Args:
        t_init (int): Transition beginning
        t_end (int): Transition end
        initvalue (int, optional): Initial value. Defaults to 0.
        endvalue (int, optional): End value. Defaults to 1.

    Returns:
        function: function that performs the linear transition
    """
    if type == 'linear':
        a = np.polyfit([t_init,t_end],[initvalue,endvalue],1)        
        f = lambda t: np.poly1d(a)(t)
        out = lambda t: initvalue*expit(10*(t_init-t)) + f(t)*(expit(10*(t-t_init)) - expit(10*(t-t_end))) + endvalue*expit(10*(t-t_end)) 

    elif type == 'quadratic':
        # convex increase or concave decrease
        if (not concavity and endvalue>=initvalue) or (concavity and initvalue>endvalue):
            t_aux = 2*t_init - t_end 
            v_aux = endvalue
            
        # concave increase or convex decrease
        elif (concavity and endvalue>=initvalue) or (not concavity and initvalue>endvalue):
            t_aux = 2*t_end - t_init 
            v_aux = initvalue
            
        a = np.polyfit([t_init,t_end,t_aux],[initvalue,endvalue,v_aux],2)
        f = lambda t: np.poly1d(a)(t)
        
        out = lambda t: initvalue*expit(10*(t_init-t)) + f(t)*(expit(10*(t-t_init)) - expit(10*(t-t_end))) + endvalue*expit(10*(t-t_end)) 
    
    elif type == 'sigmoidal':
        out = lambda t:  initvalue + (endvalue-initvalue)*(expit((t-(t_init+t_end)/2)*gw/(t_end-t_init))) 
    
    
    return out



# Saturation function
def saturation(upperlimit):
    """ Function that builds binary functions which indicates when the sum of the arguments are bigger than the saturation function.

    Args:
        upperlimit (function or cv19function builder arg): Upper limit function

    Returns:
        saturationfunction (function): binary function with time multiple arguments that returns 1 when the arguments addition function
            Args:
                t (float): time value
                *args: multiple arguments which are added and then compared with the saturation value at time t
            Returns:
                int: 
                    0 when the functions given are smaller than the saturation function 
                    1 when they are bigger or equal
    """
    gw = 20 # gain weight    
    upperlimit = build(upperlimit)

    def aux(t,*functions):
        """binary function with time multiple arguments that returns 1 when the arguments addition function

        Args:
            t (float): time value
            *args: multiple arguments which are added and then compared with the saturation value at time t
        
        Returns:
        int: 
            0 when the functions given are smaller than the saturation function 
            1 when they are bigger or equal
        """
        f = build_add(*functions)
        return(expit(gw*(f(t)-upperlimit(t))))
    return aux
    

def data_function(data,future):
    """Creates a function that returns the data during its length

    Args:
        data (list| np.array|pd.Series|pd.DataFrame ): Data 
        future ([type]): [description]
    """
    if isinstance(data, pd.DataFrame):
        data = list(data.iloc[0])
    elif isinstance(data, pd.Series):
        data = list(data)

    auxf = build(future)
        
    def aux(t):
        if t < len(data):
            return data[int(t)]
        else:
            return auxf(t)
    
    return aux


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

def interpolate_data(data):
    """Interpolates the data in order to have a daily array of data

    Args:
        data ([type]): [description]

    Returns:
        [type]: [description]
    """

    # Work in progress
    aux = data
    return aux