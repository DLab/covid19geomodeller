#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

"""
# ------------------------------------------------- #   
#                                                   #
#          Mobility Networks Construction           #
#                                                   #
# ------------------------------------------------- #

To Do:
* Build functions for:
    * Small World Flux Matrix

"""
def rnd_flux_matrix(population,fraction = 0.1):
    """Generate a random flux matrix that moves the specified fraction of the population
    Based in the dirichlet distribution in order to mantain
    
    ref: https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.dirichlet.html#numpy.random.Generator.dirichlet

    Args:
        population (array): array with the population of each town that makes part of this meta-population system. 
        fraction (float, array, optional): Population fraction that travels per day. This can also be an array that specifies a different fraction per population. Defaults to 0.1.

    Returns:
        np.array: Flux matrix with 0 diagonals 
    """
    size = len(population)
    if not isinstance(fraction,(list,np.ndarray)):
        fraction = np.ones(size)*fraction
    
    aux = []
    for i in range(size):
        aux.append(np.insert(np.random.default_rng().dirichlet(np.ones(size-1),size=1)*population[i]*fraction[i],i,0))
    return np.array(aux).astype(int)

def to_symmetric_function(inputmatrix):
    """Convert a flux matrix into a daily symmetric flux function, 
    in order to conserve the population throughout the day, avoiding long-term mass migrations

    Args:
        inputmatrix (np.array): Base flux matrix
    
    Returns:
        function: Time symmetric flux function
    """
    def Phi(t):
        if t%1<0.5:
            return inputmatrix
        else:
            return inputmatrix.transpose()
    return Phi


def rnd_flux(population,fraction=0.1):
    aux = rnd_flux_matrix(population,fraction)
    return to_symmetric_function(aux)