#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import json

"""
# ------------------------------------------------- #   
#                                                   #
#          Mobility Networks Construction           #
#                                                   #
# ------------------------------------------------- #

To Do:
* Build functions for:
    * Weighted Small World Network -> Flux Matrix
    * Weighted prefferential attachment network -> Flux matrix

"""

def gravity_model(populations, distances, beta):
    mobility_matrix = np.zeros((len(populations), len(populations)))

    for i in range(len(populations)):
        for j in range(len(populations)):
            if i != j:
                mobility_matrix[i, j] = (populations[i] * populations[j]) / (distances[i, j] ** beta)

    return mobility_matrix

def radiation_model(populations, distances, s):
    mobility_matrix = np.zeros((len(populations), len(populations)))

    for i in range(len(populations)):
        for j in range(len(populations)):
            if i != j:
                numerator = populations[i] * populations[j]
                denominator = (populations[i] + s[i, j]) * (populations[i] + populations[j] + s[i, j])
                mobility_matrix[i, j] = numerator / denominator

    return mobility_matrix


def random_mobility_model(populations, distances, fraction):
    mobility_matrix = np.zeros((len(populations), len(populations)))

    for i in range(len(populations)):
        for j in range(len(populations)):
            if i != j:
                mobility_matrix[i, j] = np.random.uniform(0, populations[i] * fraction)

        # Normalize row to ensure the total moved population is a fraction of the region's population
        mobility_matrix[i, :] /= np.sum(mobility_matrix[i, :])
        mobility_matrix[i, :] *= populations[i] * fraction

    return mobility_matrix


def create_mobility_matrix(populations, distances, model, **kwargs):
    if model == 'gravity':
        return gravity_model(populations, distances, kwargs.get('beta', 2))
    elif model == 'radiation':
        return radiation_model(populations, distances, kwargs.get('s', np.zeros_like(distances)))
    elif model == 'random':
        return random_mobility_model(populations, distances, kwargs.get('fraction', 0.1))
    else:
        raise ValueError(f"Unsupported model: {model}")
    

"""
Flux matrix generation
"""
def rnd_flux_matrix(population,fraction = 0.1,seed=None):
    """Generate a random flux matrix that moves the specified fraction of the population.
    This method uses the dirichlet distribution in order to distribute the population mantaining the total.
        
    ref: https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.dirichlet.html#numpy.random.Generator.dirichlet

    Args:
        population (array): array with the population of each town that makes part of this meta-population system. 
        fraction (float, array, optional): Population fraction that travels per day. This can also be an array that specifies a different 
        fraction per population. Defaults to 0.1.

    Returns:
        np.array: Flux matrix with 0 diagonals 
    """
    size = len(population)
    if not isinstance(fraction,(list,np.ndarray)):
        fraction = np.ones(size)*fraction
    
    aux = []
    for i in range(size):
        aux.append(np.insert(np.random.default_rng(np.random.default_rng(seed=seed)).dirichlet(np.ones(size-1),size=1)*population[i]*fraction[i],i,0))
    return np.array(aux).astype(int)


def to_symmetric_function(inputmatrix,transposed=False):
    """Convert a flux matrix into a daily symmetric flux function, 
    in order to conserve the population throughout the day, avoiding long-term mass migrations

    Args:
        inputmatrix (np.array): Base flux matrix
        transposed (bool, optional): Returns the transposed matrix
    
    Returns:
        function: Time symmetric flux function
    """
    inputmatrix_T = inputmatrix.transpose()
    def Phi(t):        
        if t%1<0.5:
            return inputmatrix
        else:
            return inputmatrix_T
    if transposed:
        def Phi_T(t):            
            if t%1<0.5:
                return inputmatrix_T
            else:
                return inputmatrix
        return Phi, Phi_T
    else:
        return Phi



def rnd_flux_symmetric(population,fraction=0.1, transposed = False,seed=None):
    return to_symmetric_function(rnd_flux_matrix(population,fraction,seed),transposed)

def rnd_flux(population,fraction=0.1,transposed=False):
    return rnd_flux_symmetric(population,fraction,transposed)


def export_mobility(mobfunction,t=None,path=None):
    aux = {}
    if not t:
        t = np.arange(0,2,0.5)
    for i in t:
        aux[i] = mobfunction(i).tolist()
    if path:
        json.dump(aux,path)
    return json.dumps(aux)

# 
def import_mobility(file):
    # Import file
    # Build Mobility Matrix
    return 

def create_mobility(inputmatrix,symmetric=True):
    # Create mobility matrix from static matrix, json    
    pass

def mobility_to_tensor(mobfunction,t_end):
    return np.array([mobfunction(i/2)for i in range(2*t_end+1)])

def mobility_transposed(matrix):
    return np.array([matrix[i].transpose() for i in range(len(matrix))])