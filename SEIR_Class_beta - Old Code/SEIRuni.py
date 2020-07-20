#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from class_SEIR import SEIR
from  SEIRrefiner import SEIRrefiner 
import pandas as pd
import numpy as np
from timeit import default_timer as timer
import Single_dist_ref_SEIR as SDSEIR
import csv
import datetime
import argparse

"""
This code calculates the optimal parameters for all the comunas using the SEIR1 method

"""

if __name__ == '__main__':
    # Get the options from the command line
    np.random.seed()
    print("SEIR Sym")
    parser = argparse.ArgumentParser(description='Dlab COVID SEIR Symulation ')
    parser.add_argument('-s', "--state", dest="state", required = False, help="State for simulation")


    args=parser.parse_args()
    if args.state:
        tstate = args.state
    else:
        tstate = None
    # Import cutlist
    cutlist = []
    cutlistpath = "../Data/cutlist.csv"
    cutlist = pd.read_csv(cutlistpath, header = None,dtype=str)

    print(tstate)
    # Input parameters:
    mov = 0.2 # movility level during quarantine 
    qp = 0 # Quarantine Period
    simdata = pd.DataFrame()
    simparameters = pd.DataFrame(index=['beta','sigma','gamma','mu','err'])
    initdate = pd.DataFrame()

    # Refine parameters per CUT    
    n = 0
    for index, row in cutlist.iterrows():    
        state = str(row[0])[0:2]
        comuna = str(row[0])
        if tstate==None or tstate == state: 
            print('Refining Comuna '+comuna)
            result = SDSEIR.ref_pygmo(state,comuna,mov,qp)
            if result:            
                simdata[comuna] = result['sim']['I'] #Check Initial date
                pd.DataFrame(result['sim']['I']).to_csv(comuna+'_I.csv')
                simparameters[comuna] = np.append(result['params'],result['err']) # Parameters
                initdate[comuna] = [(result['init_date'] + datetime.timedelta(days=43830)).strftime("%d-%b-%Y")] #Check Initial date
                print('error: '+str(result['err']))
    print('Saving csv data')
    path = '../Data/DatosConsolidados/'
    if tstate == None:
        simdata.to_csv(path+'uni_sim.csv')
        simparameters.to_csv(path+'uni_parameters.csv')
        initdate.to_csv(path+'uni_initdate.csv')
    else:
        simdata.to_csv(path+tstate+'_uni_sim.csv')
        simparameters.to_csv(path+tstate+'_uni_parameters.csv')
        initdate.to_csv(path+tstate+'_uni_initdate.csv')

    # Falta procesar los datos de salida para tener 1 por dia

    
    # {'Ir':Ir,'tr':tr, 'params':xopt, 'err':fopt,'sim':sim, 'init_date':b_date}
    