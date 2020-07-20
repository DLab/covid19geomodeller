#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:39:49 2020

@author: pmaldona
"""
import Single_dist_ref_SEIR as SDSEIR
import multiprocessing
from joblib import Parallel, delayed
import numpy as np

tsim=400
mov=0.2
rep=5


def rep_fun_eta(mu_i,rep):
    out=[] 
    for i in range(int(rep)):  
        results = SDSEIR.ref_sim_national_mu(mov=mov,qp=0,tsim = tsim,tci=None,movfunct='sawtooth',mu=mu_i)
        params=results['params']
        params=np.append(params,mu_i)
        params=np.append(params,results['err'])
        params=np.append(params,results['rerr'])
        print(params)
        out.append(params)
    return(out)   

num_cores = multiprocessing.cpu_count()

#params=Parallel(n_jobs=num_cores, verbose=50)(delayed(ref_test.refinepso_all)(Ir,tr,swarmsize=200,maxiter=50,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,1],Q_r=[0,1],obj_func='IN')for i in range(int(rep)))
params=Parallel(n_jobs=num_cores, verbose=50)(delayed(rep_fun_eta)(i,rep)for i in np.linspace(0.5,5.5,num_cores))
params=np.vstack( params )k