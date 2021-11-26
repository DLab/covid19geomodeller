#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIR Model param finder
Implementation of a Metropolis-Hasting model
"""
import class_SEIR as S
import numpy as np
import pandas
from numpy import linalg as LA

def mesh(model,b_r,s_r,g_r,m_r,Ir,tr,Npoints,steps,r0,err):
    params = []
    for i in range(Npoints):
        print(i)
        beta_i=np.random.uniform(min(),max(b_r))
        sigma_i=np.random.uniform(min(s_r),max(s_r))
        gamma_i=np.random.uniform(min(g_r),max(g_r))
        mu_i=np.random.uniform(min(g_r),max(g_r))
        print([beta_i,sigma_i,gamma_i,mu_i])
        params.append(met_hast(model,Ir,tr,beta_i,sigma_i,gamma_i,mu_i,r0,steps,err))
    return params


def met_hast(model,Ir,tr,beta_i,sigma_i,gamma_i,mu_i,r0,steps,err):
    x=S.SEIR(model.P,model.eta,model.eta,model.S[0,:],model.E[:,0],model.I[:,0],model.R[:,0],beta_i,gamma_i,sigma_i)
    x.integr_RK4(model.t0,model.T,model.h,beta_i,sigma_i,gamma_i,mu_i,False,True)
    e_0=objective_funct(Ir,tr,x.I,x.t,2)
    e_o=e_0
    params = [[beta_i,sigma_i,gamma_i,mu_i,e_o]]
    for i in range(steps):
        print(i)
        [b_p,s_p,g_p,m_p]=transition_model(x.beta,x.sigma,x.gamma,x.mu,r0,model)
        x.integr_RK4(model.t0,model.T,model.h,beta_i,sigma_i,gamma_i,mu_i,False,True)

        x_new=S.SEIR(   model.t0,model.T,model.h,b_p,s_p,g_p,m_p,False,True)
        e_n=objective_funct(Ir,tr,x_new["I"],x_new["t"],2)
        # Acceptance
        if(e_n/e_o<1):
            x=x_new
            e_o = e_n
            params.append([b_p,s_p,g_p,m_p,e_n])
        if(e_n<err):
            break

    return(params)

# objective function to minimize for any cases
def objective_funct(Ir,tr,I,t,l):
    idx=np.searchsorted(t,tr)
    return LA.norm(Ir-I[:,idx],l)


#The tranistion model defines how to move from current to new parameters
def transition_model(beta,sigma,gamma,mu,r0,model):
    rb=r0*(max(b_r)-min(b_r))
    rg=r0*(max(g_r)-min(g_r))
    rs=r0*(max(s_r)-min(s_r))
    rm=r0*(max(g_r)-min(g_r))
    b_p = np.random.normal(beta,rb)
    s_p = np.random.normal(sigma,rs)
    g_p = np.random.normal(gamma,rg)
    m_p = np.random.normal(mu,rm)

    while b_p > max(b_r) or b_p < min(b_r):
        b_p = np.random.normal(beta,rb)
    print(b_p)
    while s_p > max(s_r) or s_p < min(s_r):
        s_p = np.random.normal(sigma,rs)
    print(s_p)
    while g_p > max(g_r) or g_p < min(g_r):
        g_p = np.random.normal(gamma,rg)
    print(g_p)
    while m_p > max(g_r) or m_p < min(g_r):
        m_p = np.random.normal(beta,rb)
    print(m_p)
    return [b_p,s_p,g_p,m_p]


#if __name__ == "__main__":
    # init model
    #model = SEIR(self,P,eta,alpha,S0,E0,I0,R0,Ir,tr,h,beta_r,gamma_r,sigma_r,mu_r)
    # Le sacaria vairables de inicializacion al modelo SeIR
    # no se usan beta, gamma, sigma, mu,y h
    # I_r se usa en todos los modelos, lo abstraeria un nivel
    #   def __init__(self,P,eta,alpha,S0,E0,I0,R0,Ir,tr,h,beta_r,gamma_r,sigma_r,mu_r):
#           n = 5 # Number of initial parameters
    # create a mesh of initial parameters


    # execute met_hast for each parameter

    # Output comparison
