#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.special import expit

"""
# ------------------------------------------------- #   
#                                                   #
#        Time dependent function constructors       #
#                                                   #
# ------------------------------------------------- #
"""

def functionAddition(functarray):
    """
    Creates the function that results adding the functions present in the function array
    Receives a list where each element is a function that receives 1 argument
    """
    def f(t):
        aux = 0
        for i in functarray:    
            aux += i(t)
        return aux
    return f


# Custom values through time
def Events(values,days):
    """
    Event creator function. Create a time dependent function that returns the values setted in the 
    values vector for the periods specified in the days vector. 
    Input:
    * values: list with the values for the different intervals
    * days: list of lists of len == 2 with the values for the function on that v.
    
    """

    # for one interval
    if type(not days[0]) == list and len(days == 2):
        days = [[days[0],days[1]]]
    
    for i in range(len(values)):
        if type(values[i]) == int or type(values[i]) == float:
            aux = values[i]
            values[i] = lambda t: aux
    # define default function
    def f(t):
        return 0

    functions = [f]
    for i in range(len(values)):
        def auxf(t,i=i):
            return values[i](t)*(expit(10*(t-days[i][0])) - expit(10*(t-days[i][1]))) 
        functions.append(auxf)

    f = functionAddition(functions)
    
    return f


    """
    days = [[3,10],[15,25],[8,17]]
    values = [1,5,3]
    t = np.array(np.arange(0,30,0.1))
    plt.plot(t,sumfunct(t))
    plt.show()
    """



def square(min_val=0,max_val=1,period=14,init=0,end=1000,off_val=0,phase='min',duty=0.5):
    # phase en fraccion del periodo
    def f(t): 
        return signal.square(t,duty)
    
    if phase == 'min':
        phi = np.pi*(2*init/period - 1.1)
    elif phase == 'max':
        phi = 2*np.pi*init/period
    else:
        phi = 2*np.pi*(init/period - phase)
    def aux(t):    
        return (expit(10*(t-init)) - expit(10*(t-end)))*((max_val-min_val)/2*(f(2*np.pi / period * t - phi))+(max_val+min_val)/2) + \
        (1-(expit(10*(t-init)) - expit(10*(t-end))))*off_val

    return aux


def sine(min_val=0,max_val=1,period=14,init=0,end=1000,off_val=0,phase='min'):
    # phase en fraccion del periodo
    def f(t): 
        return np.cos(t)
    
    if phase == 'min':
        phi = np.pi*(2*init/period - 1)
    elif phase == 'max':
        phi = 2*np.pi*init/period
    else:
        phi = 2*np.pi*(init/period - phase)
    def aux(t):    
        return (expit(10*(t-init)) - expit(10*(t-end)))*((max_val-min_val)/2*(f(2*np.pi / period * t - phi))+(max_val+min_val)/2) + \
        (1-(expit(10*(t-init)) - expit(10*(t-end))))*off_val

    return aux

def sawtooth(min_val=0,max_val=1,period=14,init=0,end=1000,off_val=0,phase='min',width=1):
    # phase en fraccion del periodo
    def f(t): 
        return signal.sawtooth(t,width)
    
    if phase == 'min':
        phi = np.pi*(2*init/period - 1.1)
    elif phase == 'max':
        phi = 2*np.pi*init/period
    else:
        phi = 2*np.pi*(init/period - phase)
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
        return maxvalue*(expit((t-t0-4)*8/(t1-t0)) - expit(df*(t-t2)))      
    return f

"""
Testing:

#def increasingfunct
t = np.array(np.arange(0,50,0.1)) 
plt.plot(t,aux(t))
plt.show() 

"""