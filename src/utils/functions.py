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
    * Crear una opción para el build donde si recibe un arreglo con diccionarios creadores de funciones, los sume con functionAddition. 
    * A fitdata agregar la opción de trabajar con objetos de nuestras clases de datos. 
"""
def build(input):
    # crear un iterador que recorra el input y cree la función a partir de un comando exec:
    # Acepta diccionarios o strings con forma de diccionario

    if type(input)==str:
        input_dict = json.loads(input)
    elif type(input)==dict:
        input_dic = input.copy()
    else:
        print('Constant value function: '+str(input))
        def out(t):
            return input
        setattr(locals()['out'],'constructor',str(input))            
        return out
    try:
        print("Executing "+input_dict['function'])
    except:
        raise SyntaxError("No function defined")
    aux = 'out=' + input_dict['function']+"("
    del input_dict['function']
    for key, value in input_dict.items():
        aux +=key+"="+str(value)+',' 
    aux +=')'
    print(aux)
    ldict={}
    exec(aux,globals(),ldict)
    out = ldict['out']
    setattr(locals()['out'],'constructor',str(input))
    return locals()['out']

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

def fitdata(time,values,degree=4,tsat=-1,endvalue = -10):
    # Aproximate Hospitals capacity:
    datamodel = np.poly1d(np.polyfit(time, values, degree))
    tsat = time[tsat]
    endvalue = np.mean(values[endvalue:])
    f_out=lambda t: datamodel(t)*(1-expit(t-tsat)) + expit(t-tsat)*endvalue
    return f_out


# Custom values through time
def Events(values,days,functions = []):
    """
    Event creator function. Create a time dependent function that returns the values setted in the 
    values vector for the periods specified in the days vector. 
    Input:
    * values: list with the values for the different intervals
    * days: list of lists of len == 2 with the values for the function on that v.
    
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

    out = functionAddition(aux_f)    
    return out


    """
    days = [[3,10],[15,25],[8,17]]
    values = [1,5,3]
    t = np.array(np.arange(0,30,0.1))
    plt.plot(t,sumfunct(t))
    plt.show()
    """



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

def sawtooth(min_val=0,max_val=1,period=14,init=0,end=1000,off_val=0,initphase='min',width=1):
    # phase en fraccion del periodo
    def f(t): 
        return signal.sawtooth(t,width)
    
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