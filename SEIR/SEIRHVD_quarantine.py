#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

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