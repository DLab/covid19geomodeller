#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import json

"""
# ------------------------------------------------- #   
#                                                   #
#          Artificial Mobility Networks             #
#                                                   #
# ------------------------------------------------- #

Mobility models used for simulating the movement of people between regions:

Gravity Model:
The gravity model, inspired by Newton's law of gravitation, is a widely used spatial interaction model in various fields, such as transportation, human migration, and trade. The model assumes that the interaction between two regions (e.g., the number of people moving between them) is directly proportional to the product of their populations (or other attributes, like economic size) and inversely proportional to some power of the distance between them.

Assumptions:

The interaction between two regions is positively correlated with their population sizes.
The interaction between two regions is negatively correlated with the distance between them.
The model assumes that the system is in equilibrium, meaning that the total number of people moving from one region to another remains constant over time.
Use Cases:

Estimating traffic flows between cities.
Modeling migration patterns between countries or regions.
Predicting trade flows between countries.
Radiation Model:
The radiation model is an alternative to the gravity model for estimating spatial interactions. It was proposed to overcome some of the limitations of the gravity model, mainly the need for calibration of parameters (alpha and beta). The radiation model is derived from a probabilistic framework and does not require any free parameters. It is based on the concept of intervening opportunities, which means that the interaction between two regions is determined by the number of opportunities (e.g., jobs) that are available in other regions located between them.

Assumptions:

People move from one region to another in search of opportunities (e.g., jobs).
The likelihood of moving to a specific region is determined by the opportunities available in all the regions located between the origin and destination regions.
Use Cases:

Estimating commuting patterns in urban areas.
Modeling human mobility patterns at different spatial scales (e.g., city, country, or global levels).
Comparing Gravity and Radiation Models:
Both models have their strengths and weaknesses. The gravity model is more flexible due to its parameters, allowing it to be calibrated to better fit specific datasets. However, this also makes it more dependent on the availability of good calibration data. The radiation model does not require calibration, making it more straightforward to apply in situations where calibration data are scarce or nonexistent. Generally, the choice between the two models depends on the specific problem and the data available.


# TODO: 
* Check the symmetry of the mobility when using the daily_mobility_pattern
* Check that the population is being conserved when using the dynamic mobility functions
* Add file i/o functionality to save and load mobility matrices
"""
   
# --------------------------------- #
#           Matrix Creation         #
# --------------------------------- # 
def create_mobility_matrix(populations, distances=None, model='gravity', **kwargs):
    """Creates mobility matrix based on the specified model and parameters.

    Args:
        populations (np.array): Array of populations for each region.
        distances (np.array, optional): Matrix of distances between regions. Defaults to None.
        model (str, optional): Mobility model to be used ('gravity', 'radiation', or 'random'). Defaults to 'gravity'.
        **kwargs (dict): Additional parameters for the specific mobility models.
    Raises:
        ValueError: Error when an unsupported model is specified.

    Returns:
        np.array: Static mobility matrix.
    """
    if distances is None:
        distances = create_random_distances_matrix(len(populations), 10, 100)

    if model == 'gravity':
        return gravity_model(populations, distances, kwargs.get('beta', 2))
    elif model == 'radiation':
        return radiation_model(populations, distances, kwargs.get('s', np.zeros_like(distances)))
    elif model == 'random':
        return random_mobility_model(populations, distances, kwargs.get('fraction', 0.1))
    else:
        raise ValueError(f"Unsupported model: {model}")


def create_dynamic_mobility(mobility_model, dynamic_pattern, populations, distances=None, **kwargs):
    """Create a dynamic mobility matrix function based on the specified mobility model and dynamic pattern.

    Args:
        mobility_model (str): Mobility model to be used ('gravity', 'radiation', or 'random')._
        dynamic_pattern (str): Dynamic pattern to be used ('sinusoidal', 'rush_hour', or 'piecewise_linear').
        populations (np.array): Array of populations for each region.
        distances (np.array, optional): Matrix of distances between regions. Defaults to None, which generates a random distance matrix.
        **kwargs (dict): Additional parameters for the specific mobility models and dynamic patterns.
    Raises:
        ValueError: Error when an unsupported dynamic pattern is specified.

    Returns:
        function: Dynamic mobility matrix function.
    """
    if distances is None:
        distances = create_random_distances_matrix(len(populations), 10, 100)

    static_mobility = create_mobility_matrix(mobility_model, populations, distances, **kwargs)

    if dynamic_pattern == 'sinusoidal':
        return daily_mobility_pattern(static_mobility, **kwargs)
    elif dynamic_pattern == 'rush_hour':
        return rush_hour_mobility_pattern(static_mobility, **kwargs)
    elif dynamic_pattern == 'piecewise_linear':
        return piecewise_linear_mobility_pattern(static_mobility, **kwargs)
    else:
        raise ValueError(f"Unknown dynamic pattern: {dynamic_pattern}")


# --------------------------------- #
#           Mobility Models         #
# --------------------------------- # 
def random_mobility_model(population, fraction=0.5, seed=None,**kwargs):
    """Generate a random flux matrix that moves the specified fraction of the population.
    This method uses the dirichlet distribution in order to distribute the population maintaining the total.

    Args:
        population (array): array with the population of each town that makes part of this meta-population system.
        fraction (float, array, optional): Population fraction that travels per day. This can also be an array that specifies a different fraction per population. Defaults to 0.5.

    Returns:
        np.array: Flux matrix with 0 diagonals
    """
    size = len(population)
    if not isinstance(fraction, (list, np.ndarray)):
        fraction = np.ones(size) * fraction

    rng = np.random.default_rng(seed=seed)
    
    aux = []
    for i in range(size):
        aux.append(np.insert(rng.dirichlet(np.ones(size - 1), size=1) * population[i] * fraction[i], i, 0))
    return np.array(aux).astype(int)

def gravity_model(populations, distances, alpha=1, beta=1, fraction=0.5):
    """Calculate the gravity model mobility matrix.
    The gravity model, inspired by Newton's law of gravitation, is a widely used spatial interaction model in various fields, such as transportation, human migration, and trade. The model assumes that the interaction between two regions (e.g., the number of people moving between them) is directly proportional to the product of their populations (or other attributes, like economic size) and inversely proportional to some power of the distance between them.
    
    Args:
        populations (np.array): Array of populations for each region.
        distances (np.array): Matrix of distances between regions.
        alpha (int, optional): Parameter to control the effect of population sizes. Defaults to 1.
        beta (int, optional):  Parameter to control the effect of distances.. Defaults to 1.
        fraction (float, optional): Fraction of the population that travels per day. Defaults to 0.5.

    Returns:
        np.array: Mobility matrix based on the gravity model.
    """
    num_regions = len(populations)
    mobility_matrix = np.zeros((num_regions, num_regions))

    for i in range(num_regions):
        for j in range(num_regions):
            if i != j:
                mobility_matrix[i, j] = (populations[i] ** alpha) * (populations[j] ** alpha) / (distances[i, j] ** beta)
    row_sums = mobility_matrix.sum(axis=1, keepdims=True)
    mobility_matrix = mobility_matrix / row_sums * populations[:, np.newaxis] * fraction
    
    return mobility_matrix.astype(int)

def radiation_model(populations, distances, fraction=0.5):
    """Calculate the radiation model mobility matrix.

    Args:
        populations (np.array): Array of populations for each region.
        distances (np.array): Matrix of distances between regions.
        fraction (float, optional): Fraction of the population that travels per day. Defaults to 0.5.

    Returns:
        np.array: Mobility matrix based on the radiation model.
    """
    num_regions = len(populations)
    mobility_matrix = np.zeros((num_regions, num_regions))

    for i in range(num_regions):
        for j in range(num_regions):
            if i != j:
                s_ij = populations[i] * populations[j]
                m_ij = (s_ij * distances[i, j]) / ((populations[i] + s_ij) * (populations[j] + s_ij))
                mobility_matrix[i, j] = m_ij

    # Normalize rows to prevent the total outgoing flux from exceeding the population
    row_sums = mobility_matrix.sum(axis=1, keepdims=True)
    mobility_matrix = mobility_matrix / row_sums * populations[:, np.newaxis] * fraction
    return mobility_matrix.astype(int)




# --------------------------------------- #
#         Dynamic Mobility Patterns       #
# --------------------------------------- # 
def daily_mobility_pattern(mobility_matrix, amplitude=0.5, period=1, phase_shift=0):
    """Create a time-varying mobility matrix function using a sinusoidal pattern.

    Args:
        mobility_matrix (np.array): Static mobility matrix.
        amplitude (float, optional): Amplitude of the sinusoidal pattern. Defaults to 0.5.
        period (int, optional): Period of the sinusoidal pattern. Defaults to 1 day.
        phase_shift (int, optional): Phase shift of the sinusoidal pattern. Defaults to 0.
    
    Returns:
        function: Time-varying mobility matrix function.
    """
    def time_varying_mobility(t):
        scaling_factor = 1 + amplitude * np.sin(2 * np.pi * (t / period - phase_shift))
        return mobility_matrix * scaling_factor

    return time_varying_mobility

def rush_hour_mobility_pattern(mobility_matrix, peak_amplitude=0.5, off_peak_amplitude=-0.5, peak_start=7/24, peak_end=9/24, evening_start=17/24, evening_end=19/24):
    """Create a time-varying mobility matrix function using a rush hour pattern.

    Args:
        mobility_matrix (np.array): Static mobility matrix.
        peak_amplitude (float, optional): Amplitude during peak hours. Defaults to 0.5.
        off_peak_amplitude (float, optional): Amplitude during off-peak hours. Defaults to -0.5.
        peak_start (float, optional): Start time of morning peak hours. Defaults to 7/24.
        peak_end (float, optional): End time of morning peak hours. Defaults to 9/24.
        evening_start (float, optional):  Start time of evening peak hours. Defaults to 17/24.
        evening_end (float, optional): End time of evening peak hours. Defaults to 19/24.
        
    Returns:
        function: Time-varying mobility matrix function.
    """
    def time_varying_mobility(t):
        if peak_start <= t < peak_end:
            scaling_factor = 1 + off_peak_amplitude + (peak_amplitude - off_peak_amplitude) * (t - peak_start) / (peak_end - peak_start)
        elif peak_end <= t < evening_start:
            scaling_factor = 1 + off_peak_amplitude
        elif evening_start <= t < evening_end:
            scaling_factor = 1 + off_peak_amplitude + (peak_amplitude - off_peak_amplitude) * (t - evening_start) / (evening_end - evening_start)
        else:
            scaling_factor = 1 + off_peak_amplitude
        return mobility_matrix * scaling_factor

    return time_varying_mobility


def piecewise_linear_mobility_pattern(mobility_matrix, intervals, slopes):
    """Create a time-varying mobility matrix function using a piecewise linear pattern.

    Args:
        mobility_matrix (np.array): Static mobility matrix
        intervals (list of tuples): List of time intervals represented as (start, end) tuples.
        slopes (list of float): List of slopes for each time interval.
        
    Returns:
        function: Time-varying mobility matrix function.
    """
    def time_varying_mobility(t):
        for i, (start, end) in enumerate(intervals):
            if start <= t < end:
                scaling_factor = 1 + slopes[i] * (t - start) / (end - start)
                break
        else:
            scaling_factor = 1

        return mobility_matrix * scaling_factor

    return time_varying_mobility


# --------------------------------- #
#           Distance Matrix         #
# --------------------------------- # 
def create_random_distances_matrix(size, min_distance=10,max_distance=100, seed=None):
    """
    Create a random symmetric distance matrix.

    Args:
        size (int): The number of nodes (regions) in the distance matrix.
        max_distance (float): The maximum distance between any two nodes (regions).
        seed (int, optional): Seed for the random number generator.

    Returns:
        np.array: A symmetric distance matrix.
    """
    rng = np.random.default_rng(seed=seed)
    distances = rng.uniform(min_distance, max_distance, size=(size, size))
    distances = (distances + distances.T) / 2
    np.fill_diagonal(distances, 0)
    return distances

# --------------------------------- #
#              File I/O             #
# --------------------------------- #
def export_mobility_dynamic(mobfunction,t=None,path=None):
    aux = {}
    if not t:
        t = np.arange(0,2,0.5)
    for i in t:
        aux[i] = mobfunction(i).tolist()
    if path:
        json.dump(aux,path)
    return json.dumps(aux)

def export_mobility_matrix(mobility_matrix,file_path):
    # Export mobility matrix to file
    return

# 
def import_mobility(file):
    # Import file
    # Build Mobility Matrix
    return 

