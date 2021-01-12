#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.special import expit

"""
# ------------------------------------------------- #   
#                                                   #
#       SEIRHDV Quarantine Scenarios Definition     #
#                                                   #
# ------------------------------------------------- #

Quarantine array:
np.array([tsim,max_mob,rem_mob,qp,iqt,fqt,movfunct])

* tsim: Tiempo de simulación
* max_mob: Movilidad máxima durante tiempo sin cuarentena
* rem_mob: Movilidad remanente durante tiempo de cuarentena
* qp: Periodo de Cuarentena para cuarentenas alternantes - qp días con cuarentena, luego qp días sin cuarentena
* iqt: Día de inicio de cuarentena (desde el inicio de la simulación)
* fqt: Día de fin de cuarentena (desde el inicio de la simulación)
* movfunct: Función de movilidad
    * 0: Cuarentena total durante el período comprendido entre iqt y fqt
    * 1: Cuarentena alternante tipo onda Cuadrada con período qp
    * 2: Cuarnetena tipo diente de cierra con período qp


"""
def functionAdd(funct):
    def f(t):
        aux = 0
        for i in funct:    
            aux += i(t)
        return aux
    return f

class Quarantine():
    """
    Quarantine Object
        input: 
            rem_mob,max_mob=0.85,qp=0,iqt=0,fqt=1000,movfunct = 'once'
        output:
            Quarantine.plot()
            Quarantine.alpha()
        Usage example:
            alpha = Quarantine(rem_mob,max_mob=0.85,qp=0,iqt=0,fqt=1000,movfunct = 'once').alpha(t)
    """
    def __init__(self,rem_mob,max_mob=0.85,qp=0,iqt=0,fqt=1000,movfunct = 'once'):
        self.rem_mob = rem_mob
        self.max_mob = max_mob
        self.qp = qp
        self.iqt = iqt
        self.fqt = fqt
        self.movfunct = movfunct

        self.alpha = self.alphafunct()


    def alphafunct(self):
        """    
        # max_mob: Movilidad sin cuarentena
        # rem_mob: Movilidad con cuarentena
        # qp: Periodo cuarentena dinamica 
        #          - qp >0 periodo Qdinamica 
        #          - qp = 0 sin qdinamica
        # iqt: Initial quarantine time. Tiempo inicial antes de cuarentena dinamica
        #          - iqt>0 inicia con cuarentena total hasta iqt
        #          - iqt<0 sin cuarentena hasta iqt
        # fqt: Final quarantine time. Duracion tiempo cuarentena 
        # movfunct: Tipo de cuarentena dinamica desde iqt
        #          - once: una vez durante qp dias 
        #          - total: total desde iqt
        #          - sawtooth: diente de cierra
        #          - square: onda cuadrada


        # FQT no funciona en las ondas periódicas
        """
        def alpha(t):             
            if 'square' in self.movfunct or self.movfunct == 1:
               def f(t): 
                   return signal.square(t)
               if t<abs(self.iqt):
                   if self.iqt>0:
                       return(self.rem_mob)
                   else:
                       return(self.max_mob)
               else:
                   if self.qp == 0:
                       return(self.max_mob)
                   elif t<self.fqt:
                       return((self.max_mob-self.rem_mob)/2*(f(np.pi / self.qp * t - np.pi))+(self.max_mob+self.rem_mob)/2)
                   else:
                       return(self.max_mob)   


            elif 'once' in self.movfunct or self.movfunct == 0:        
                if t<self.iqt:
                    return(self.max_mob)
                elif t>self.fqt:
                    return(self.max_mob)
                else:
                    return(self.rem_mob)


            elif 'sawtooth' in self.movfunct or self.movfunct == 2:
               def f(t): 
                   return signal.sawtooth(t)
               if t<abs(self.iqt):
                   if self.iqt>0:
                       return(self.rem_mob)
                   else:
                       return(self.max_mob)
               else:
                   if self.qp == 0:
                       return(self.max_mob)
                   elif t<self.fqt:
                       return((self.max_mob-self.rem_mob)/2*(f(np.pi / self.qp * t - np.pi))+(self.max_mob+self.rem_mob)/2)
                   else:
                       return(self.max_mob)   
                     
        return(alpha)

    def plot(self,endtime=100):
        t = list(range(endtime))
        mobility = [self.alpha(i) for i in t]
        plt.plot(list(t),mobility,label='Mobility function')
        #plt.legend(loc=0)
        plt.show()



def alphafunct(rem_mob,max_mob=0.85,qp=0,iqt=0,fqt=1000,movfunct = 'once'):
    """    
    # max_mob: Movilidad sin cuarentena
    # rem_mob: Movilidad con cuarentena
    # qp: Periodo cuarentena dinamica 
    #          - qp >0 periodo Qdinamica 
    #          - qp = 0 sin qdinamica
    # iqt: Initial quarantine time. Tiempo inicial antes de cuarentena dinamica
    #          - iqt>0 inicia con cuarentena total hasta iqt
    #          - iqt<0 sin cuarentena hasta iqt
    # fqt: Final quarantine time. Duracion tiempo cuarentena 
    # movfunct: Tipo de cuarentena dinamica desde iqt
    #          - once: una vez durante qp dias 
    #          - total: total desde iqt
    #          - sawtooth: diente de cierra
    #          - square: onda cuadrada
    """
    def alpha(t):             
        if 'square' in movfunct or movfunct == 1:
            def f(t): 
                return signal.square(t)
            if t<abs(iqt):
                if iqt>0:
                    return(rem_mob)
                else:
                    return(max_mob)
            else:
                if qp == 0:
                    return(max_mob)
                elif t<fqt:
                    return((max_mob-rem_mob)/2*(f(np.pi / qp * t - np.pi))+(max_mob+rem_mob)/2)
                else:
                    return(max_mob)   


        elif 'once' in movfunct or movfunct == 0:        
            if t<iqt:
                return(max_mob)
            elif t>fqt:
                return(max_mob)
            else:
                return(rem_mob)


        elif 'sawtooth' in movfunct or movfunct == 2:
            def f(t): 
                return signal.sawtooth(t)
            if t<abs(iqt):
                if iqt>0:
                    return(rem_mob)
                else:
                    return(max_mob)
            else:
                if qp == 0:
                    return(max_mob)
                elif t<fqt:
                    return((max_mob-rem_mob)/2*(f(np.pi / qp * t - np.pi))+(max_mob+rem_mob)/2)
                else:
                    return(max_mob)   
                    
    return(alpha)


# Increasing function
def SeroPrevDynamics(t0,t1,t2,dailyincrease = 1,form='sig',df = 10):
    """
    Sero Prevalence Dynamics Function generator
    """
    if form == 'line':
        a = np.polyfit([t0,t1],[0,dailyincrease],1)
        def f(t):
            return(np.poly1d(a)(t))

        def chi(t):
            return f(t)*(expit(10*(t-t0)) - expit(10*(t-t1))) + dailyincrease*(expit(10*(t-t1)) - expit(10*(t-t2)))
    elif form == 'quadratic':
        a = np.polyfit([t0,t1,2*t0-t1],[0,dailyincrease,dailyincrease],2)
        def f(t):
            return(np.poly1d(a)(t))

        def chi(t):
            return f(t)*(expit(10*(t-t0)) - expit(10*(t-t1))) + dailyincrease*(expit(10*(t-t1)) - expit(10*(t-t2)))            
    elif form == 'sig'or form == 'sigmoid':
        def chi(t):
            return dailyincrease*(expit((t-t0-4)*8/(t1-t0)) - expit(df*(t-t2)))
      
    return chi


"""
def functionSum(a,b):
    def aux(t):
        return a(t)+b(t)
    return aux
"""


# square function
def Exams(examrate,period,duty=50):
    """
    Examination campain 
    """
    def psi(t):
        return (examrate*signal.square(2*np.pi*t/period,duty) + examrate)/2
  
    return psi


# Custom values through time
def Events(values,days):
    """
    Event creator function. Create a time dependent function that returns the values setted in the 
    values vector for the periods specified in the days vector. 
    Input:
    * values: list with the values for the different intervals
    * days: list of lists of len == 2 with the values for the function on that v.
    By default the value returned is 0 unless is stated otherwise
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

    f = functionAdd(functions)
    
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


"""
Testing:

#def increasingfunct
t = np.array(np.arange(0,50,0.1)) 
plt.plot(t,aux(t))
plt.show() 

"""