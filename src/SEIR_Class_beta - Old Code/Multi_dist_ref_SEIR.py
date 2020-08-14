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
import multiprocessing
from joblib import Parallel, delayed
import json
import requests
import datetime as dt
from scipy import signal

def refine_multi(state,comunas):
        
    #Endpoint encuesta SETRA, no funcional aún
    #    endpoint="http://192.168.2.223:5004/v1/get_movilidad_por_comuna"
    #    r = requests.get(endpoint) #, params = {"w":"774508"})
    #    mydict = r.json()
    #    OD=pd.DataFrame(mydict)
        
    #Endpoint encuensta ENE    
    #    endpoint="http://192.168.2.223:5004/v2/get_movilidad?month=1&year=2019"
    #    r = requests.get(endpoint) #, params = {"w":"774508"})
    #    mydict = r.json()
    #    info=pd.DataFrame(mydict)
        
    # Data csv encuenstra Sectra (funcional)
    path="../Data/UrbanCenters/"
    OD=pd.read_csv(path+"Movilidad_diaria_comuna.csv",index_col=0)
    P=OD.filter(comunas) 
    P=P.loc[comunas].to_numpy()
    
    endpoint="http://192.168.2.220:8080/covid19/findComunaByIdState?idState="+state+"&&comuna="
    r = requests.get(endpoint) #, params = {"w":"774508"})
    mydict = r.json()
    info=pd.DataFrame(mydict)
    
    comunas=info[info.description.isin(comunas)]
    
    init_dates=[]
    for i in comunas.cut:
        endpoint="http://192.168.2.220:8080/covid19/getDatosMinSalSummary?state="+state+"&comuna="+i
        r = requests.get(endpoint) #, params = {"w":"774508"})
        mydict = r.json()
        data=pd.DataFrame(mydict)
        data=data[data.data != 0]
        init_dates.append( dt.datetime.strptime(data.labels.iloc[0], '%d/%m'))
    
    init_date=max(init_dates).strftime('%d/%m')
    
    Ic=[]
    for i in comunas.cut:
        endpoint="http://192.168.2.220:8080/covid19/getDatosMinSalSummary?state="+state+"&comuna="+i
        r = requests.get(endpoint) #, params = {"w":"774508"})
        mydict = r.json()
        data=pd.DataFrame(mydict)
        init_i=data[data.labels==init_date].index
        if len(init_i)!=0 :
            init_i=data[data.labels==init_date].index[0]
            if len(Ic)==0:
                Ic=data[['labels','data']][init_i:].set_index('labels')
                Ic.rename(columns={'data': i}, inplace=True)
            else:
                Ic=Ic.join(data[['labels','data']][init_i:-1].set_index('labels'), on='labels')
                Ic.rename(columns={'data': i}, inplace=True)
        else:
            init_i=0
            if len(Ic)==0:
                Ic=data[['labels','data']][init_i:].set_index('labels')
                Ic.rename(columns={'data': i}, inplace=True)
            else:
                Ic=Ic.join(data[['labels','data']][init_i:-1].set_index('labels'), on='labels')
                Ic.rename(columns={'data': i}, inplace=True)
    
    Ic=Ic.fillna(0)
    
    for i in (Ic.columns):
        for j in range(max(Ic.count())-1):
            if Ic[i].iloc[j+1]<Ic[i].iloc[j]:
                Ic[i].iloc[j+1]=Ic[i].iloc[j]
    
    
    tr=np.zeros(int(max(Ic.count())))
    for i in range(1,int(max(Ic.count())),1):
        diff=dt.datetime.strptime(Ic.index[i], '%d/%m')-dt.datetime.strptime(Ic.index[i-1], '%d/%m')
        tr[i]=diff.days+tr[i-1]
    
    Cn=np.zeros(Ic.shape)
    for i in range(Ic.shape[1]):
        Cn[0][i]=Ic.iloc[0][i]
        for j in range(1,Ic.shape[0],1):
            Cn[j][i]=Ic.iloc[j][i]-Ic.iloc[j-1][i]
    
    Ir=np.zeros(Ic.shape)
    for i in range(Ic.shape[1]):
        Ir[np.where(tr-14<=0)[0],i]=Ic.iloc[np.where(tr-14<=0)[0],i]
        for j in np.where(tr-14>0)[0]:
            ind=np.where(tr-tr[j]+14<0)[0]            
            Ir[j,i]=Ic.iloc[j,i]-sum(Cn[0:ind[-1],i])
     
        
    So=comunas.numPopulation.to_numpy()
    Ro=np.zeros(P.shape[0])
    Ir = Ic.transpose().to_numpy()
    Io = Ir[:,0]
    Eo = 2*Io

    S0 = So -Eo -Io
    E0 = Eo
    I0 = Io
    R0 = Ro

    h=0.05

    #b_r=[0.01,3.5] #0.2
    b_r=0.2 #0.2
    s_r=[0.16,0.25] #  5 dias
    g_r=[0.03,0.14] # 14 dias
    mu_r=[0.5,3.5] #2

    # Strategy functions
    # In first instance the 
    def alpha(t):
        return(np.zeros([Ic.shape[1],Ic.shape[1]]))

    def eta(t):
        return(np.ones(Ic.shape[1]))
      
    #Multiprocressing refine method via PSO
    
    rep=5 #number of repetition by mu parameter
    params=[]
    
    # Mu optimization for #rep repetitions 
    def rep_fun_eta(mu,rep):
        out=[] 
        for i in range(int(rep)):  
            #ref_test.refinepso_all(Ir,tr,swarmsize=100,maxiter=50,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,10],Q_r=[0,10],obj_func='IR')
            ref_test=SEIRrefiner(P,eta,alpha,S0,E0,I0,R0,min(tr),max(tr),h,b_r,s_r,g_r,mu)
            ref_test.refinepso_eta(Ir,tr,swarmsize=200,maxiter=100,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,10],obj_func='IN')
            out.append(ref_test.paramsPSO)
        return(out)   
    
    num_cores = multiprocessing.cpu_count()
    
    #params=Parallel(n_jobs=num_cores, verbose=50)(delayed(ref_test.refinepso_all)(Ir,tr,swarmsize=200,maxiter=50,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,1],Q_r=[0,1],obj_func='IN')for i in range(int(rep)))
    params=Parallel(n_jobs=num_cores, verbose=50)(delayed(rep_fun_eta)(i,rep)for i in np.linspace(min(mu_r),max(mu_r),num_cores))
    params=np.vstack( params )
    
    ind=np.where(params[:,-1]==min(params[:,-1]))
    alpha_test=False
    if len(ind[0])>0: 
        for i in range(len(ind[0])):
            cand=np.random.choice(ind[0])
            if np.all(params[cand,4:Ic.shape[1]+4]!=0):
                params_e=params[cand,:]
                alpha_test=True
                print(params_e)
                break
            elif(len(ind[0])>1):
                ind[0]=np.delete(ind[0],cand)
            else:
                params_e=params[cand,:]
                break
    if not alpha_test:

        f_params=params_e[0:Ic.shape[1]+4]
        print(f_params)
        f_params=np.append(f_params,np.zeros(Ic.shape[1]))
        print(f_params)
        f_params=np.append(f_params,params_e[Ic.shape[1]+5:Ic.shape[1]+6])
        
    else:
        # Update eta function
        def eta(t):
            return(params_e[4:Ic.shape[1]+4])
        
        def rep_fun_alpha(mu,rep):
            out=[] 
            for i in range(int(rep)):  
                #ref_test.refinepso_all(Ir,tr,swarmsize=100,maxiter=50,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,10],Q_r=[0,10],obj_func='IR')
                ref_test=SEIRrefiner(P,eta,alpha,S0,E0,I0,R0,min(tr),max(tr),h,b_r,s_r,g_r,mu)
                ref_test.refinepso_alpha(Ir,tr,swarmsize=200,maxiter=100,omega=0.5, phip=0.5, phig=0.5,Q_r=[0,1],obj_func='IN')
                out.append(ref_test.paramsPSO)
            return(out)
    
        params=Parallel(n_jobs=num_cores, verbose=50)(delayed(rep_fun_alpha)(i,rep)for i in np.linspace(min(mu_r),max(mu_r),24))
        params=np.vstack( params )
    
    
        ind=np.where(params[:,-1]==min(params[:,-1]))
        print("ind inicial")
        print(ind)
        if len(ind[0])>0: 
            for i in range(len(ind[0])):
                cand=np.random.choice(ind[0])
                if np.any(params[cand,4:Ic.shape[1]+4]!=0):
                    params_a=params[cand,:]
                    alpha_test=True
                    print("ind encontrado")
                    print(ind)
                    break
                elif(len(ind[0])>1):
                    ind[0]=np.delete(ind[0],cand)
                else:
                    params_a=params[cand,:]
                    break
        print(params_a)
        print(params_a[0:4])
        f_params=params_a[0:4]
        print(f_params,params_e[4:Ic.shape[1]+4])
        f_params=np.append(f_params,params_e[4:Ic.shape[1]+4])
        print(f_params,params_a[4:Ic.shape[1]+6])
        f_params=np.append(f_params,params_a[4:Ic.shape[1]+6])
        
    return({'params':f_params,'init_date':init_date})
    
def simulate_multi(state,comunas,params,qp=0,mov=0.2,tsim=300,tci=None,movfunct='sawtooth'):
    # params = beta sigma gamma mu
    # Movility function shape
    if movfunct == 'sawtooth':
        def f(t): 
            return signal.sawtooth(t)
    elif movfunct == 'square':
        def f(t):
            return signal.square(t)
    
    if qp <=0:
        def e_fun(t):
            return(1)
    
    elif tci is None:        
        def e_fun(t):
            return((1-mov)/2*(f(np.pi / qp * t - np.pi))+(1+mov)/2)
    else:
        def e_fun(t):
                if t<tci:
                    return(1)
                else:
                    return((1-mov)/2*(f(np.pi / qp * t - np.pi))+(1+mov)/2)
    

    #Endpoint encuesta SECTR, no funcional aún
    #    endpoint="http://192.168.2.223:5004/v1/get_movilidad_por_comuna"
    #    r = requests.get(endpoint) #, params = {"w":"774508"})
    #    mydict = r.json()
    #    OD=pd.DataFrame(mydict)
        
    #Endpoint encuensta ENE    
    #    endpoint="http://192.168.2.223:5004/v2/get_movilidad?month=1&year=2019"
    #    r = requests.get(endpoint) #, params = {"w":"774508"})
    #    mydict = r.json()
    #    info=pd.DataFrame(mydict)
        
    # Data csv encuenstra Sectra (funcional)
    path="../Data/"
    OD=pd.read_csv(path+"Movilidad_diaria_comuna.csv",index_col=0)
    P=OD.filter(comunas) 
    P=P.loc[comunas].to_numpy()
    
    endpoint="http://192.168.2.220:8080/covid19/findComunaByIdState?idState="+state+"&&comuna="
    r = requests.get(endpoint) #, params = {"w":"774508"})
    mydict = r.json()
    info=pd.DataFrame(mydict)
    
    comunas=info[info.description.isin(comunas)]
    
    init_dates=[]
    for i in comunas.cut:
        endpoint="http://192.168.2.220:8080/covid19/getDatosMinSalSummary?state="+state+"&comuna="+i
        r = requests.get(endpoint) #, params = {"w":"774508"})
        mydict = r.json()
        data=pd.DataFrame(mydict)
        data=data[data.data != 0]
        init_dates.append( dt.datetime.strptime(data.labels.iloc[0], '%d/%m'))
    
    init_date=max(init_dates).strftime('%d/%m')
    
    Ic=[]
    for i in comunas.cut:
        endpoint="http://192.168.2.220:8080/covid19/getDatosMinSalSummary?state="+state+"&comuna="+i
        r = requests.get(endpoint) #, params = {"w":"774508"})
        mydict = r.json()
        data=pd.DataFrame(mydict)
        init_i=data[data.labels==init_date].index
        if len(init_i)!=0 :
            init_i=data[data.labels==init_date].index[0]
            if len(Ic)==0:
                Ic=data[['labels','data']][init_i:].set_index('labels')
                Ic.rename(columns={'data': i}, inplace=True)
            else:
                Ic=Ic.join(data[['labels','data']][init_i:-1].set_index('labels'), on='labels')
                Ic.rename(columns={'data': i}, inplace=True)
        else:
            init_i=0
            if len(Ic)==0:
                Ic=data[['labels','data']][init_i:].set_index('labels')
                Ic.rename(columns={'data': i}, inplace=True)
            else:
                Ic=Ic.join(data[['labels','data']][init_i:-1].set_index('labels'), on='labels')
                Ic.rename(columns={'data': i}, inplace=True)
    
    Ic=Ic.fillna(0)
    
    for i in (Ic.columns):
        for j in range(max(Ic.count())-1):
            if Ic[i].iloc[j+1]<Ic[i].iloc[j]:
                Ic[i].iloc[j+1]=Ic[i].iloc[j]
    
    
    tr=np.zeros(int(max(Ic.count())))
    for i in range(1,int(max(Ic.count())),1):
        diff=dt.datetime.strptime(Ic.index[i], '%d/%m')-dt.datetime.strptime(Ic.index[i-1], '%d/%m')
        tr[i]=diff.days+tr[i-1]
    
    Cn=np.zeros(Ic.shape)
    for i in range(Ic.shape[1]):
        Cn[0][i]=Ic.iloc[0][i]
        for j in range(1,Ic.shape[0],1):
            Cn[j][i]=Ic.iloc[j][i]-Ic.iloc[j-1][i]
    
    Ir=np.zeros(Ic.shape)
    for i in range(Ic.shape[1]):
        Ir[np.where(tr-14<=0)[0],i]=Ic.iloc[np.where(tr-14<=0)[0],i]
        for j in np.where(tr-14>0)[0]:
            ind=np.where(tr-tr[j]+14<0)[0]            
            Ir[j,i]=Ic.iloc[j,i]-sum(Cn[0:ind[-1],i])
     
        
    So=comunas.numPopulation.to_numpy()
    Ro=np.array([0,0,0,0,0])
    Ir = Ic.transpose().to_numpy()
    Io = Ir[:,0]
    Eo = 2*Io

    S0 = So -Eo -Io
    E0 = Eo
    I0 = Io
    R0 = Ro

    h=0.01
    
    ia=3+P.shape[0]
    def alpha(t):
        alpha=np.zeros([P.shape[0],P.shape[0]])
        for j in range(P.shape[0]):
             for k in range(P.shape[0]):
                 alpha[j][k]=e_fun(t)*params[j+ia]*params[k+ia]
        return(alpha)

    def eta(t):
        eta=np.ones(P.shape[0])
        for j in range(P.shape[0]):
            eta[j]=e_fun(t)*params[3+j]
        return(eta)

    
    sim = SEIR(P,eta,alpha,S0,E0,I0,R0,params[0],params[1],params[2],params[3])
    sim.integr(min(tr),tsim,h,True)

    return({'S':sim.S,'E':sim.E,'I':sim.I,'R':sim.R,'tout':sim.t,'init_date':init_date})
    # Cada elemento es una matriz con las comunas en las columnas y el tiempo en las filas. El orden de las comunas
    # Es el mismo que entra como argumento. 
    
