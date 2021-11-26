#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from class_SEIR import SEIR
from  SEIRrefiner import SEIRrefiner 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from timeit import default_timer as timer

if __name__ == '__main__':
    start = timer()
    # Get the options from the command line
    np.random.seed()
    print("SEIR Sym")
    parser = argparse.ArgumentParser(description='Dlab COVID SEIR Symulation ')
    parser.add_argument('-o', "--outputfile", dest="outputfile", required = True, help="Output file to store data")
    parser.add_argument('-pr', "--paramrange", dest="prange", required = False, help="Parameter Range")

    args=parser.parse_args()
    outfile = args.outputfile
    #print("Output File: %s" ,% args.outputfile)

    method = "MH"

    """
    #---------------------#
    #     Import Data     #
    #---------------------#
    """
    if os.name == 'nt':
        path = "D:\\samue\\Dropbox\\AFES Datascience\\Ciencia y Vida\\Data\\"

    elif os.name == 'posix':
        path = "../Data/"

    So = pd.read_excel(path+"poblacion_Inicial_S_stgo.xlsx", header=None).to_numpy()    
    Eo = pd.read_excel(path+"poblacion_Inicial_E_stgo.xlsx", header=None).to_numpy()    
    Io = pd.read_excel(path+"poblacion_Inicial_I_stgo.xlsx", header=None).to_numpy()    
    Ro = pd.read_excel(path+"poblacion_Inicial_R_stgo.xlsx", header=None).to_numpy()
    P = pd.read_excel(path+"connectivity_stgo2.xlsx", header=None).to_numpy()
    Ir = pd.read_excel(path+"Simulacion-400dias-I.xlsx", header=None).to_numpy()
    S_su1 = pd.read_excel(path+"Simulacion-400dias-S.xlsx", header=None).to_numpy()
    E_su1 = pd.read_excel(path+"Simulacion-400dias-E.xlsx", header=None).to_numpy()
    I_su1 = pd.read_excel(path+"Simulacion-400dias-I.xlsx", header=None).to_numpy()
    R_su1 = pd.read_excel(path+"Simulacion-400dias-R.xlsx", header=None).to_numpy()


    # Create model variables
    S0 = So[:,0]
    E0 = Eo[:,0]
    I0 = Io[:,0]
    R0 = Ro[:,0]

    #Init variables
    tr=np.arange(Ir.shape[1])
    h=0.1

    # Parameter range 
    # If not parameter range file do this:
    b_r=[0.1,0.3] #0.2
    s_r=[0.05,0.15] #0.1
    g_r=[0.05,0.15] #0.1
    mu_r=[1,3] #2


    # Strategy functions
    def alpha(t):
        return(np.ones([34,34])-np.eye(34))

    def eta(t):
        return(np.ones(34))

    # Object init

    """
    #---------------------------------------------
    #   Find paramters using different Methos    #
    #---------------------------------------------
    """
    # ---------------------------- #
    # Create param refiner object  #
    # ---------------------------- #
    ref_test=SEIRrefiner(P,eta,alpha,S0,E0,I0,R0,min(tr),max(tr),0.1,b_r,g_r,s_r,mu_r)

    # ---------------------------- #
    #     Metropolis-Hastings      #
    # ---------------------------- # 
    if "MH" in method:
        tr=np.arange(Ir.shape[1])   
        ref_test.refineMH(Ir,tr,0.1,1,6400,1)
        print("listeilor")
        print(ref_test.paramsMH[0])
        # Save parameters 
        print("Saving data")    
        np.savetxt(outfile, ref_test.paramsMH[0], delimiter=",")

    # ------------------------------ #
    #        Random Walk             #
    # ------------------------------ # 
    if "RW" in method:
        tr=np.arange(Ir.shape[1])   
        ref_test.refineRW(Ir,tr,0.1,2,2,1)
        print("listeilor")
        print(ref_test.paramsRW)
        # Save parameters 
        print("Saving data")    
        np.savetxt(outfile, ref_test.paramsRW, delimiter=",")

    
    # ---------------------------------#
    #   Particle Swarm Optimization    #
    # -------------------------------- # 
    # beta mu y gamma
    if "PSO" in method:
        tr=np.arange(Ir.shape[1])  
        ref_test.refinepso(Ir,tr,swarmsize=2,maxiter=2,omega=0.5, phip=0.5, phig=0.5)
        print(ref_test.paramsPSO)
        print("listeilor")
        print(ref_test.paramsPSO)
        # Save parameters 
        print("Saving data")    
        np.savetxt(outfile, ref_test.paramsPSO, delimiter=",")        


    stop = timer()
    time = stop - start
    with open("time.txt", "w") as f:
        f.write(str(time))
    print(str(time))
    # SEIR Object test
    # test = SEIR(P,eta,alpha,S0,E0,I0,R0,np.mean(b_r),np.mean(g_r),np.mean(s_r),np.mean(mu_r))
    # test.integr_RK4(min(tr),max(tr),0.1,True)


    # # mesh test

    # # parms=mesh(test,Ir,tr,5,5,0.025,20)
    # idx=np.searchsorted(test.S,tr)

    # test.S[1,idx]-S_su1[1,:]




    # idx=np.searchsorted(test.t,tr)

    # S_dif=S_su1-test.S[:,idx]
    # E_dif=E_su1-test.E[:,idx]
    # I_dif=I_su1-test.I[:,idx]
    # R_dif=R_su1-test.R[:,idx]

    # np.amax(S_dif)
    # np.amax(E_dif)
    # np.amax(I_dif)
    # np.amax(R_dif)

    # plt.figure()
    # plt.plot(test.t[idx],S_dif[1,:],label='Susceptible')
    # plt.plot(test.t[idx],E_dif[1,:],label='Exposed')
    # plt.plot(test.t[idx],I_dif[1,:],label='Infected simulation')
    # plt.plot(test.t[idx],R_dif[1,:],label='Removed')
    # plt.xlabel('Days')
    # plt.ylabel('Population')
    # plt.title('COVID-19 Model')
    # plt.legend(loc=0)
    # plt.show()

    # plt.figure()
    # plt.plot(test.t[idx],test.S[1,idx],label='Susceptible')
    # plt.plot(test.t[idx],test.E[1,idx],label='Exposed')
    # plt.plot(test.t[idx],test.I[1,idx],label='Infected simulation')
    # plt.plot(test.t[idx],test.R[1,idx],label='Removed')
    # plt.xlabel('Days')
    # plt.ylabel('Population')
    # plt.title('COVID-19 Model')
    # plt.legend(loc=0)
    # plt.show()


    # plt.figure()
    # plt.plot(test.t[idx],S_su1[1,:],label='Susceptible')
    # plt.plot(test.t[idx],E_su1[1,:],label='Exposed')
    # plt.plot(test.t[idx],I_su1[1,:],label='Infected simulation')
    # plt.plot(test.t[idx],R_su1[1,:],label='Removed')
    # plt.xlabel('Days')
    # plt.ylabel('Population')
    # plt.title('COVID-19 Model')
    # plt.legend(loc=0)
    # plt.show()


    #run mesh

    #run mesh



    # Report results
