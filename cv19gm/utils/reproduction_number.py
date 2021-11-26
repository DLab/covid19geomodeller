# -*- coding: utf-8 -*-
"""
Fundación Ciencia y Vida
dlab : Computational Biology Lab
Questions to : sropert@dlab.cl, fcastillo@dlab.cl
"""
import pandas as pd
# rpy2 imports
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects
pandas2ri.activate()

# R library
global eps
eps = importr("EpiEstim")



def R_effective(smoothed, window):

    df = pd.DataFrame({"dates": smoothed.index.tolist(), "cases": smoothed.values})
    rdf = pandas2ri.py2ri(df)
    r = robjects.r
    window = window  # always 5, test
    start = r.seq(2, df.values.shape[0]-window)
    end = r.seq( 2+ window,df.values.shape[0])
    results = eps.estimate_R(rdf[1], method="parametric_si", config = eps.make_config(t_start = start,
                               t_end = end,
                               mean_si = 5,
                               std_si = 2))
    results2 = dict(results.items())
    rhat = pandas2ri.ri2py(results2["R"])
    rhat.index = smoothed.index[2+window-1:]
    return rhat

def Interpolate(Data, time):
    pass

def R_basic(data):
    pass