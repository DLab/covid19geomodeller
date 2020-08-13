#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal

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

class SEIRHVD_quarantine():
    #------------------------------------------------- #
    #              Definir Escenarios                  #
    #------------------------------------------------- #
    #tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct
    def defaultescenarios(self):
        self.inputarray=np.array([
                [self.tsim,0.85,0.6,0,self.May15,500,0],
                [self.tsim,0.85,0.65,0,self.May15,500,0],
                [self.tsim,0.85,0.7,0,self.May15,500,0]])        
        self.numescenarios = len(self.inputarray)

    def addquarantine(self,tsim=None,max_mov=None,rem_mov=None,qp=None,iqt=None,fqt=None,movfunct=None):        
        if tsim:
            self.inputarray.append([tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct])
        self.numescenarios= len(self.inputarray)
        return()

    #def addquarantinevector(self,tsim=None,max_mov=None,rem_mov=None,qp=None,iqt=None,fqt=None,movfunct=None):        
    #    if tsim:
    #        self.inputarray.append([tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct])
    #    self.numescenarios= len(self.inputarray)
    #    return()        

    # traspasarlo a lenguaje humao
    def showscenarios(self):        
        print(self.inputarray)
        return()


    def setscenario(self,tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct):
        self.tsim = tsim
        self.max_mov = max_mov
        self.rem_mov = rem_mov
        self.qp = qp
        self.iqt = iqt
        self.fqt = fqt
        if movfunct == 0:
            self.movfunct = 'once'
        elif movfunct == 1:
            self.movfunct = 'square'
        elif movfunct == 2:
            self.movfunct = 'sawtooth'
        else:
            self.movfunct = 'once'
        
        self.alpha = self.alphafunct(self.max_mov,self.rem_mov,self.qp,self.iqt,self.fqt,self.movfunct)
        return() 
        


    def alphafunct(self,max_mov,rem_mov,qp,iqt=0,fqt=300,movfunct = 'once'):
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
            if 'square' in movfunct:
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


            elif 'once' in movfunct:        
                if t<iqt:
                    return(max_mov)
                elif t>fqt:
                    return(max_mov)
                else:
                    return(rem_mov)


            elif 'sawtooth' in movfunct:
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



    def get_conditions(self):
        # This function will return a text explaining the different simulation scenarios currently placed as inputs            
        return
