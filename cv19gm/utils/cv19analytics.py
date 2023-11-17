from datetime import timedelta
import numpy as np

"""Analytics module

"""

def peak(model):
    """
    Perform simulation analytics after running it.
    It calculates peaks, prevalence, and will include R(t). 
    """
    #Cálculo de la fecha del Peak  
    model.peakindex = np.where(model.I==max(model.I))[0][0]
    model.peak = max(model.I)
    model.peak_t = model.t[model.peakindex]
    if model.initdate:
        model.dates = [model.initdate+timedelta(int(model.t[i])) for i in range(len(model.t))]
        model.peak_date = model.initdate+timedelta(days=int(round(model.peak_t)))
    else:
        model.dates = [None for i in range(len(model.t))]
        model.peak_date = None            

    # Prevalence: 
    model.prevalence_total = model.I_ac/model.population
    model.prevalence_susc = [model.I_ac[i]/(model.S[i]+model.E[i]+model.I[i]+model.R[i]) for i in range(len(model.I_ac))]
    model.prevalence_det = [model.pI_det*model.I_ac[i]/(model.S[i]+model.E[i]+model.I[i]+model.R[i]) for i in range(len(model.I_ac))]                         
    return
        
        
def calculateindicators(model):
    model.R_ef
    model.SHFR

    # SeroPrevalence Calculation

    # Errors (if real data)

    # Active infected
    print('work in progress')