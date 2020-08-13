#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

"""
# ------------------------------------------------- #   
#                                                   #
#       SEIRHDV Quarantine Scenarios Definition     #
#                                                   #
# ------------------------------------------------- #

Quarantine array:
np.array([tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct])

* tsim: Tiempo de simulación
* max_mov: Movilidad máxima durante tiempo sin cuarentena
* rem_mov: Movilidad remanente durante tiempo de cuarentena
* qp: Periodo de Cuarentena para cuarentenas alternantes - qp días con cuarentena, luego qp días sin cuarentena
* iqt: Día de inicio de cuarentena (desde el inicio de la simulación)
* fqt: Día de fin de cuarentena (desde el inicio de la simulación)
* movfunct: Función de movilidad
    * 0: Cuarentena total durante el período comprendido entre iqt y fqt
    * 1: Cuarentena alternante tipo onda Cuadrada con período qp
    * 2: Cuarnetena tipo diente de cierra con período qp


"""

class Quarantine():
    def __init__(self,rem_mov,max_mov=0.85,qp=0,iqt=0,fqt=1000,movfunct = 'once'):
        self.rem_mov = rem_mov
        self.max_mov = max_mov
        self.qp = qp
        self.iqt = iqt
        self.fqt = fqt
        self.movfunct = movfunct

        self.alpha = self.alphafunct


    def alphafunct(self):
        """    
        # max_mov: Movilidad sin cuarentena
        # rem_mov: Movilidad con cuarentena
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
        def alpha(self,t):             
            if 'square' in self.movfunct or self.movfunct == 1:
               def f(t): 
                   return signal.square(t)
               if t<abs(self.iqt):
                   if self.iqt>0:
                       return(self.rem_mov)
                   else:
                       return(self.max_mov)
               else:
                   if self.qp == 0:
                       return(self.max_mov)
                   elif t<self.fqt:
                       return((self.max_mov-self.rem_mov)/2*(f(np.pi / qp * t - np.pi))+(self.max_mov+self.rem_mov)/2)
                   else:
                       return(self.max_mov)   


            elif 'once' in self.movfunct or self.movfunct == 0:        
                if t<self.iqt:
                    return(self.max_mov)
                elif t>self.qt:
                    return(self.max_mov)
                else:
                    return(self.rem_mov)


            elif 'sawtooth' in self.movfunct or self.movfunct == 2:
               def f(t): 
                   return signal.sawtooth(t)
               if t<abs(self.iqt):
                   if self.iqt>0:
                       return(self.rem_mov)
                   else:
                       return(self.max_mov)
               else:
                   if self.qp == 0:
                       return(self.max_mov)
                   elif t<self.fqt:
                       return((self.max_mov-self.rem_mov)/2*(f(np.pi / qp * t - np.pi))+(self.max_mov+self.rem_mov)/2)
                   else:
                       return(self.max_mov)   
                     
        return(alpha)

    def plot(self,endtime):
        t = range(endtime)
        mobility = [self.alpha(i) for i in t]
        plt.plot(t,mobility,legend='Mobility function')
        plt.show()



def alphafunct(rem_mov,max_mov=0.85,qp=0,iqt=0,fqt=1000,movfunct = 'once'):
    """    
    # max_mov: Movilidad sin cuarentena
    # rem_mov: Movilidad con cuarentena
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
                    return(rem_mov)
                else:
                    return(max_mov)
            else:
                if qp == 0:
                    return(max_mov)
                elif t<fqt:
                    return((max_mov-rem_mov)/2*(f(np.pi / qp * t - np.pi))+(max_mov+rem_mov)/2)
                else:
                    return(max_mov)   


        elif 'once' in movfunct or movfunct == 0:        
            if t<iqt:
                return(max_mov)
            elif t>fqt:
                return(max_mov)
            else:
                return(rem_mov)


        elif 'sawtooth' in movfunct or movfunct == 2:
            def f(t): 
                return signal.sawtooth(t)
            if t<abs(iqt):
                if iqt>0:
                    return(rem_mov)
                else:
                    return(max_mov)
            else:
                if qp == 0:
                    return(max_mov)
                elif t<fqt:
                    return((max_mov-rem_mov)/2*(f(np.pi / qp * t - np.pi))+(max_mov+rem_mov)/2)
                else:
                    return(max_mov)   
                    
    return(alpha)
