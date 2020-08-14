#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIR Model param finder
Implementation of a Metropolis-Hasting model
"""
from class_SEIR import SEIR
import numpy as np
from numpy import linalg as LA
from pyswarm import pso

class SEIRrefiner:
    """
    constructor of SEIR class elements, it's initialized when a parameter
    miminization is performed to adjust the best setting of the actual infected
    """
    def __init__(self,P,eta,alpha,S0,E0,I0,R0,t0,T,h,beta_r,sigma_r,gamma_r,mu_r):

        #init response values
        self.S0=S0
        self.E0=E0
        self.I0=I0
        self.R0=R0
        self.N=(S0+E0+I0+R0).astype("float64")
        self.beta_r=beta_r
        self.sigma_r=sigma_r
        self.gamma_r=gamma_r
        self.mu_r=mu_r
        self.P = P
        
        # Output Params
        
        # Metropolis-Hastings  
        self.paramsMH = []
        self.optimalMH = []
        self.errorMH = []

        # Random Walk
        self.paramsRW = []
        self.optimalRW = []
        self.errorRW = []

        # PSO
        self.paramsPSO = []
        self.optimalPSO = []
        self.errorPSO = []

        self.error=None
        self.SEIR=None

        #saved stragegy functions
        self.alpha=alpha
        self.eta=eta

        #definition of timegrid
        self.t0=t0
        self.T=T
        self.h=h
        self.t=np.arange(self.t0,self.T+self.h,self.h)


    """
    # ---------------------------- #
    #   Params Refine Strategies   #
    # ---------------------------- #
    """
    def refineMH(self,I_r,tr,r0,Npoints,steps,err,obj_func = 'IN'):

        if obj_func == 'IN':
            objective_function = self.objective_funct_I_Norm
        elif obj_func == 'IR':
            objective_function = self.objective_funct_I_rate            
        # find the optimal parameter using metropolis-hastings
        # desnity of phase space parameter for the optimizacion of calculeted infected vs. real infected
        mesh = self.mesh(Npoints)
        # tr=np.arange(I_r.shape[1])
        results = []
        for i in range(Npoints):            
            aux = self.met_hast(I_r,tr,mesh[i][0],mesh[i][1],mesh[i][2],mesh[i][3],r0,steps,err,objective_function)
            results.append(aux)
            
        self.paramsMH = np.array(results)
        optindex = np.where(self.paramsMH[:,4]==np.amin(self.paramsMH[:,4]))[0][0]
        self.optimalMH =self.paramsMH[optindex,:]
        self.errorMH=self.optimalMH[-1]
        # should define an exit protocol
        self.SEIR = SEIR(self.P,self.eta,self.alpha,self.S0,self.E0,self.I0,self.R0,self.optimalMH[0],self.optimalMH[1],self.optimalMH[2],self.optimalMH[3])
        return self.optimalMH

    def refineRW(self,I_r,tr,r0,Npoints,steps,err,obj_func = 'IN'):
        # Refine using Montecarlo - Markov Chain
        mesh = self.mesh(Npoints)
        
        if obj_func == 'IN':
            objective_function = self.objective_funct_I_Norm
        elif obj_func == 'IR':
            objective_function = self.objective_funct_I_rate

        results = []
        for i in range(Npoints):
            # print("Mesh point number "+str(i))
            aux = self.LocMinRW(I_r,tr,mesh[i][0],mesh[i][1],mesh[i][2],mesh[i][3],r0,steps,err,objective_function)
            results.append(aux)
            # print("Error: "+str(aux[-1]))
        self.paramsRW = np.array(results)
        optindex = np.where(self.paramsRW[:,4]==np.amin(self.paramsRW[:,4]))[0][0]
        self.optimalRW =self.paramsRW[optindex,:]
        self.errorRW=self.optimalRW[-1]
        # should define an exit protocol
        self.SEIR = SEIR(self.P,self.eta,self.alpha,self.S0,self.E0,self.I0,self.R0,self.optimalRW[0],self.optimalRW[1],self.optimalRW[2],self.optimalRW[3])
        return self.optimalRW

    def refinepso(self,Ir,tr,swarmsize=5,maxiter=25,omega=0.5, phip=0.5, phig=0.5,obj_func='IN'):
        if obj_func == 'IN':
            objective_function = self.objective_funct_I_Norm
        elif obj_func == 'IR':
            objective_function = self.objective_funct_I_rate

        self.paramsPSO = self.pso_opt(Ir,tr,objective_function,omega, phip, phig,swarmsize,maxiter)
        self.optimalRW = self.paramsPSO
        self.errorPSO = self.paramsPSO[-1]
        return self.paramsPSO 

    def refinepso_all(self,Ir,tr,swarmsize=5,maxiter=25,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,10],Q_r=[0,10],obj_func='IN'):
        if obj_func == 'IN':
            objective_function = self.objective_funct_I_Norm
        elif obj_func == 'IR':
            objective_function = self.objective_funct_I_rate

        self.paramsPSO = self.pso_opt_all(Ir,tr,objective_function,omega, phip, phig,swarmsize,maxiter,eta_r=[0,2],Q_r=[0,1])
        self.optimalRW = self.paramsPSO
        self.errorPSO = self.paramsPSO[-1]
        return self.paramsPSO 
    
    def refinepso_eta(self,Ir,tr,swarmsize=5,maxiter=25,omega=0.5, phip=0.5, phig=0.5,eta_r=[0,10],obj_func='IN'):
        if obj_func == 'IN':
            objective_function = self.objective_funct_I_Norm
        elif obj_func == 'IR':
            objective_function = self.objective_funct_I_rate

        self.paramsPSO = self.pso_opt_eta(Ir,tr,objective_function,omega, phip, phig,swarmsize,maxiter,eta_r=[0,2])
        self.optimalRW = self.paramsPSO
        self.errorPSO = self.paramsPSO[-1]
        return self.paramsPSO
    
    def refinepso_alpha(self,Ir,tr,swarmsize=5,maxiter=25,omega=0.5, phip=0.5, phig=0.5,Q_r=[0,10],obj_func='IN'):
        if obj_func == 'IN':
            objective_function = self.objective_funct_I_Norm
        elif obj_func == 'IR':
            objective_function = self.objective_funct_I_rate

        self.paramsPSO = self.pso_opt_alpha(Ir,tr,objective_function,omega, phip, phig,swarmsize,maxiter,Q_r=[0,1])
        self.optimalRW = self.paramsPSO
        self.errorPSO = self.paramsPSO[-1]
        return self.paramsPSO 
    """
    # ------------------------------- #
    #   Params Optimization Methods   #
    # ------------------------------- #
    """
    

    
    # Particle Swarm Method
    def pso_opt_all(self,Ir,tr,objective_funct,omega=0.5, phip=0.5, phig=0.5,swarmsize=5,maxiter=25,eta_r=[0,2],Q_r=[0,1]):
        # _(self,P,eta,alpha,S0,E0,I0,R0,beta,gamma,sigma,mu):
        # def integr(self,t0,T,h,E0init=False):  
        
        lb = []
        ub = []

        # Beta
        if type(self.beta_r) == list:            
            lb.append(min(self.beta_r))
            ub.append(max(self.beta_r))  
        # Sigma
        if type(self.sigma_r) == list:            
            lb.append(min(self.sigma_r))
            ub.append(max(self.sigma_r))                    
        # Gamma
        if type(self.gamma_r) == list:
            lb.append(min(self.gamma_r))
            ub.append(max(self.gamma_r))                
        # Mu
        if type(self.mu_r) == list:
            lb.append(min(self.mu_r))
            ub.append(max(self.mu_r))                
        
        for j in range(self.P.shape[0]):
            lb.append(min(eta_r))
            ub.append(max(eta_r))
        
        for j in range(self.P.shape[0]):
            lb.append(min(Q_r))
            ub.append(max(Q_r))
              
        
        def opti(x):
            i = 0
            aux =[0,0,0,0]            
            # Beta
            if type(self.beta_r) == list:
                aux[0] = x[i]
                i+=1                
            else: 
                aux[0] = self.beta_r            
            # Sigma
            if type(self.sigma_r) == list:
                aux[1] = x[i]                 
                i+=1                                
            else: 
                aux[1] = self.sigma_r            
            # Gamma
            if type(self.gamma_r) == list:
                aux[2] = x[i]             
                i+=1                                
            else: 
                aux[2] = self.gamma_r
            # Mu
            if type(self.mu_r) == list:
                aux[3] = x[i]              
                i+=1                                
            else: 
                aux[3] = self.mu_r                
            ie=0
            ie=i
            def eta(t):
                eta=np.ones(self.P.shape[0])
                for j in range(self.P.shape[0]):
                    eta[j]=x[ie+j]
                return(eta)
            self.eta=eta
            
            ia=0
            ia=ie+self.P.shape[0]

            def alpha(t):
                alpha=np.zeros([self.P.shape[0],self.P.shape[0]])
                for j in range(self.P.shape[0]):
                     for k in range(self.P.shape[0]):
                         alpha[j][k]=x[j+ia]*x[k+ia]
                return(alpha)

            self.alpha=alpha
            
            model = SEIR(self.P,eta,alpha,self.S0,self.E0,self.I0,self.R0,aux[0],aux[1],aux[2],aux[3])
            model.integr(self.t0,self.T,self.h,True)
            self.eta=eta
            self.alpha=alpha

            return(objective_funct(Ir,tr,model.I,model.t,'fro')) 

        
        xopt, fopt = pso(opti, lb, ub, minfunc=1e-8, omega=omega, phip=phip, phig=phig,swarmsize=swarmsize,maxiter=maxiter)
        
        aux = [0,0,0,0]
        i = 0
        if type(self.beta_r) == list:
            aux[0] = xopt[i]
            i+=1                
        else: 
            aux[0] = self.beta_r
        
        # Sigma
        if type(self.sigma_r) == list:
            aux[1] = xopt[i]
            i+=1                                
        else: 
            aux[1] = self.sigma_r   
        
        # Gamma
        if type(self.gamma_r) == list:
            aux[2] = xopt[i]              
            i+=1                                
        else: 
            aux[2] = self.gamma_r

        # Mu
        if type(self.mu_r) == list:
            aux[3] = xopt[i]
            i+=1                                                  
        else: 
            aux[3] = self.mu_r                
            
        for j in range(len(xopt)-i):
            aux=np.append(aux,xopt[j+i])
            
        aux = np.append(aux,fopt)
        aux = np.append(aux,fopt/LA.norm(Ir,'fro')) #Porcentual error
        
        return aux

        # Particle Swarm Method
    def pso_opt_eta(self,Ir,tr,objective_funct,omega=0.5, phip=0.5, phig=0.5,swarmsize=5,maxiter=25,eta_r=[0,2]):
        # _(self,P,eta,alpha,S0,E0,I0,R0,beta,gamma,sigma,mu):
        # def integr(self,t0,T,h,E0init=False):  
        
        lb = []
        ub = []

        # Beta
        if type(self.beta_r) == list:            
            lb.append(min(self.beta_r))
            ub.append(max(self.beta_r))  
        # Sigma
        if type(self.sigma_r) == list:            
            lb.append(min(self.sigma_r))
            ub.append(max(self.sigma_r))                    
        # Gamma
        if type(self.gamma_r) == list:
            lb.append(min(self.gamma_r))
            ub.append(max(self.gamma_r))                
        # Mu
        if type(self.mu_r) == list:
            lb.append(min(self.mu_r))
            ub.append(max(self.mu_r))                
        
        for j in range(self.P.shape[0]):
            lb.append(min(eta_r))
            ub.append(max(eta_r))
        

        def opti(x):
            i = 0
            aux =[0,0,0,0]            
            # Beta
            if type(self.beta_r) == list:
                aux[0] = x[i]
                i+=1                
            else: 
                aux[0] = self.beta_r            
            # Sigma
            if type(self.sigma_r) == list:
                aux[1] = x[i]                 
                i+=1                                
            else: 
                aux[1] = self.sigma_r            
            # Gamma
            if type(self.gamma_r) == list:
                aux[2] = x[i]             
                i+=1                                
            else: 
                aux[2] = self.gamma_r
            # Mu
            if type(self.mu_r) == list:
                aux[3] = x[i]              
                i+=1                                
            else: 
                aux[3] = self.mu_r                
            
            ie=0
            ie=i
            def eta(t):
                eta=np.ones(self.P.shape[0])
                for j in range(self.P.shape[0]):
                    eta[j]=x[ie+j]
                return(eta)
            
            model = SEIR(self.P,eta,self.alpha,self.S0,self.E0,self.I0,self.R0,aux[0],aux[1],aux[2],aux[3])
            model.integr(self.t0,self.T,self.h,True)
            self.eta=eta

            return(objective_funct(Ir,tr,model.I,model.t,'fro')) 

        # Refine Eta - Alpha is fixed 
        xopt, fopt = pso(opti, lb, ub, minfunc=1e-8, omega=omega, phip=phip, phig=phig,swarmsize=swarmsize,maxiter=maxiter)
        
        aux = [0,0,0,0]
        i = 0
        if type(self.beta_r) == list:
            aux[0] = xopt[i]
            i+=1                
        else: 
            aux[0] = self.beta_r
        
        # Sigma
        if type(self.sigma_r) == list:
            aux[1] = xopt[i]
            i+=1                                
        else: 
            aux[1] = self.sigma_r   
        
        # Gamma
        if type(self.gamma_r) == list:
            aux[2] = xopt[i]              
            i+=1                                
        else: 
            aux[2] = self.gamma_r

        # Mu
        if type(self.mu_r) == list:
            aux[3] = xopt[i]
            i+=1                                                  
        else: 
            aux[3] = self.mu_r                
            
        for j in range(len(xopt)-i):
            aux=np.append(aux,xopt[j+i])
            
        aux = np.append(aux,fopt)
        aux = np.append(aux,fopt/LA.norm(Ir,'fro')) #Porcentual error
        
        return aux

    
    
    def pso_opt_alpha(self,Ir,tr,objective_funct,omega=0.5, phip=0.5, phig=0.5,swarmsize=5,maxiter=25,Q_r=[0,1]):
        # _(self,P,eta,alpha,S0,E0,I0,R0,beta,gamma,sigma,mu):
        # def integr(self,t0,T,h,E0init=False):  
        
        lb = []
        ub = []

        # Beta
        if type(self.beta_r) == list:            
            lb.append(min(self.beta_r))
            ub.append(max(self.beta_r))  
        # Sigma
        if type(self.sigma_r) == list:            
            lb.append(min(self.sigma_r))
            ub.append(max(self.sigma_r))                    
        # Gamma
        if type(self.gamma_r) == list:
            lb.append(min(self.gamma_r))
            ub.append(max(self.gamma_r))                
        # Mu
        if type(self.mu_r) == list:
            lb.append(min(self.mu_r))
            ub.append(max(self.mu_r))                
        
        
        for j in range(self.P.shape[0]):
            lb.append(min(Q_r))
            ub.append(max(Q_r))
              
        
        def opti(x):
            i = 0
            aux =[0,0,0,0]            
            # Beta
            if type(self.beta_r) == list:
                aux[0] = x[i]
                i+=1                
            else: 
                aux[0] = self.beta_r            
            # Sigma
            if type(self.sigma_r) == list:
                aux[1] = x[i]                 
                i+=1                                
            else: 
                aux[1] = self.sigma_r            
            # Gamma
            if type(self.gamma_r) == list:
                aux[2] = x[i]             
                i+=1                                
            else: 
                aux[2] = self.gamma_r
            # Mu
            if type(self.mu_r) == list:
                aux[3] = x[i]              
                i+=1                                
            else: 
                aux[3] = self.mu_r                

            ia=i
            def alpha(t):
                alpha=np.zeros([self.P.shape[0],self.P.shape[0]])
                for j in range(self.P.shape[0]):
                     for k in range(self.P.shape[0]):
                        if j==k:
                            alpha[j][k]=0
                        else:         
                            alpha[j][k]=x[j+ia]*x[k+ia]
                return(alpha)

            
            model = SEIR(self.P,self.eta,alpha,self.S0,self.E0,self.I0,self.R0,aux[0],aux[1],aux[2],aux[3])
            model.integr(self.t0,self.T,self.h,True)
            self.alpha=alpha

            return(objective_funct(Ir,tr,model.I,model.t,'fro')) 

        
        xopt, fopt = pso(opti, lb, ub, minfunc=1e-8, omega=omega, phip=phip, phig=phig,swarmsize=swarmsize,maxiter=maxiter)
        
        aux = [0,0,0,0]
        i = 0
        if type(self.beta_r) == list:
            aux[0] = xopt[i]
            i+=1                
        else: 
            aux[0] = self.beta_r
        
        # Sigma
        if type(self.sigma_r) == list:
            aux[1] = xopt[i]
            i+=1                                
        else: 
            aux[1] = self.sigma_r   
        
        # Gamma
        if type(self.gamma_r) == list:
            aux[2] = xopt[i]              
            i+=1                                
        else: 
            aux[2] = self.gamma_r

        # Mu
        if type(self.mu_r) == list:
            aux[3] = xopt[i]
            i+=1                                                  
        else: 
            aux[3] = self.mu_r                
            
        for j in range(len(xopt)-i):
            aux=np.append(aux,xopt[j+i])
            
        aux = np.append(aux,fopt)
        aux = np.append(aux,fopt/LA.norm(Ir,'fro')) #Porcentual error
        
        return aux

    # Particle Swarm Method
    def pso_opt(self,Ir,tr,objective_funct,omega=0.5, phip=0.5, phig=0.5,swarmsize=5,maxiter=25):
        # _(self,P,eta,alpha,S0,E0,I0,R0,beta,gamma,sigma,mu):
        # def integr(self,t0,T,h,E0init=False):  
        
        lb = []
        ub = []

        # Beta
        if type(self.beta_r) == list:            
            lb.append(min(self.beta_r))
            ub.append(max(self.beta_r))  
        # Sigma
        if type(self.sigma_r) == list:            
            lb.append(min(self.sigma_r))
            ub.append(max(self.sigma_r))                    
        # Gamma
        if type(self.gamma_r) == list:
            lb.append(min(self.gamma_r))
            ub.append(max(self.gamma_r))                
        # Mu
        if type(self.mu_r) == list:
            lb.append(min(self.mu_r))
            ub.append(max(self.mu_r))                

        # print(lb)
        # print(ub)            
      
        def opti(x):
            i = 0
            aux =[0,0,0,0]            
            # Beta
            if type(self.beta_r) == list:
                aux[0] = x[i]
                i+=1                
            else: 
                aux[0] = self.beta_r            
            # Sigma
            if type(self.sigma_r) == list:
                aux[1] = x[i]                 
                i+=1                                
            else: 
                aux[1] = self.sigma_r            
            # Gamma
            if type(self.gamma_r) == list:
                aux[2] = x[i]             
                i+=1                                
            else: 
                aux[2] = self.gamma_r
            # Mu
            if type(self.mu_r) == list:
                aux[3] = x[i]              
                i+=1                                
            else: 
                aux[3] = self.mu_r                

            model = SEIR(self.P,self.eta,self.alpha,self.S0,self.E0,self.I0,self.R0,aux[0],aux[1],aux[2],aux[3])
            model.integr(self.t0,self.T,self.h,True)
            

            return(objective_funct(Ir,tr,model.I,model.t,'fro')) 

        
        xopt, fopt = pso(opti, lb, ub, minfunc=1e-8, omega=omega, phip=phip, phig=phig,swarmsize=swarmsize,maxiter=maxiter)
        
        aux = [0,0,0,0]
        i = 0
        if type(self.beta_r) == list:
            aux[0] = xopt[i]
            i+=1                
        else: 
            aux[0] = self.beta_r
        
        # Sigma
        if type(self.sigma_r) == list:
            aux[1] = xopt[i]
            i+=1                                
        else: 
            aux[1] = self.sigma_r   
        
        # Gamma
        if type(self.gamma_r) == list:
            aux[2] = xopt[i]              
            i+=1                                
        else: 
            aux[2] = self.gamma_r

        # Mu
        if type(self.mu_r) == list:
            aux[3] = xopt[i]
            i+=1                                                  
        else: 
            aux[3] = self.mu_r                

        aux = np.append(aux,fopt)
        aux = np.append(aux,fopt/LA.norm(Ir,'fro')) #Porcentual error
        
        return aux

    
    
    # Metropolis Hastings
    def met_hast(self,I_r,tr,beta_i,sigma_i,gamma_i,mu_i,r0,steps,err,objective_funct):
        #print("Build SEIR")
        x=SEIR(self.P,self.eta,self.alpha,self.S0,self.E0,self.I0,self.R0,beta_i,gamma_i,sigma_i,mu_i)
        #print("RK4")
        if x.scikitsimport:
            x.integr(self.t0,self.T,self.h,True)
        else:
            #print("RK4")
            x.integr_RK4(self.t0,self.T,self.h,True)
        e_0=objective_funct(I_r,tr,x.I,x.t,2)
        e_o=e_0
        params = [[beta_i,sigma_i,gamma_i,mu_i,e_o]]
        i=0
        #print("Met-Hast")
        while i <steps:
            [b_p,s_p,g_p,m_p]=self.transition_model(x.beta,x.sigma,x.gamma,x.mu,r0)
            x_new=SEIR(self.P,self.eta,self.alpha,self.S0,self.E0,self.I0,self.R0,b_p,s_p,g_p,m_p)
            if x_new.scikitsimport:
                x_new.integr(self.t0,self.T,self.h,True)
            else:
                x_new.integr_RK4(self.t0,self.T,self.h,True)
            e_n=objective_funct(I_r,tr,x_new.I,x_new.t,2)
     
            # Acceptance
            if(e_n/e_o<1):
                x=x_new
                e_o = e_n
                params.append([b_p,s_p,g_p,m_p,e_n])
                i+=1
                print("------------------")
                print(b_p,s_p,g_p,m_p)
                print(e_o,e_n)
                # print("time: "+str(end-start))
                # print("Step "+str(i))
            else:
                u=np.random.uniform(0,1)
                if(e_n/e_0>u):
                    x=x_new
                    e_o = e_n
                    params.append([b_p,s_p,g_p,m_p,e_n])
                    i+=1
                    print("------------------")
                    print(b_p,s_p,g_p,m_p)
                    print(e_o,e_n)
        #sleep(0.01)
        # Dejo los params historicos por si hay que debuggear esta parte
        return params


    # Local Minimum Random Walk"
    # Like a RW with bias towards the optimum
    def LocMinRW(self,I_r,tr,beta_i,sigma_i,gamma_i,mu_i,r0,steps,err,objective_funct):
        #print("Build SEIR")
        x=SEIR(self.P,self.eta,self.alpha,self.S0,self.E0,self.I0,self.R0,beta_i,gamma_i,sigma_i,mu_i)
        #print("RK4")
        if x.scikitsimport:
            x.integr(self.t0,self.T,self.h,True)
        else:
            #print("RK4")
            x.integr_RK4(self.t0,self.T,self.h,True)
        e_0=objective_funct(I_r,tr,x.I,x.t,2)
        e_o=e_0
        params = [[beta_i,sigma_i,gamma_i,mu_i,e_o]]
        i=0
        k=0
        #print("Met-Hast")
        while i <steps:
            [b_p,s_p,g_p,m_p]=self.transition_model(x.beta,x.sigma,x.gamma,x.mu,r0)
            x_new=SEIR(self.P,self.eta,self.alpha,self.S0,self.E0,self.I0,self.R0,b_p,s_p,g_p,m_p)
            if x_new.scikitsimport:
                x_new.integr(self.t0,self.T,self.h,True)
            else:
                x_new.integr_RK4(self.t0,self.T,self.h,True)
            e_n=objective_funct(I_r,tr,x_new.I,x_new.t,2)
     
            # Acceptance
            if(e_n/e_o<1):
                x=x_new
                e_o = e_n
                params.append([b_p,s_p,g_p,m_p,e_n])
                i+=1
                k=0
                print("------------------")
                print(b_p,s_p,g_p,m_p)
                print(e_o,e_n)
                # print("time: "+str(end-start))
                # print("Step "+str(i))
            if(e_n<err):
                break
            k+=1
            if k>=100:
                break
            #sleep(0.01)
        # Dejo los params historicos por si hay que debuggear esta parte
        return params[-1]
    
    """
    # ------------------------------------ #
    #   Optimization Objective Functions   #
    # ------------------------------------ #
    Arguments: Ir: Data Matrix, tr: Time vector, I: Simulated Data, t , 
    """
    # objective function to minimize Infected difference norm
    def objective_funct_I_Norm(self,Ir,tr,I,t,l):
        idx=np.searchsorted(t,tr)
        return LA.norm(Ir-I[:,idx],l)

    
    # objective function to minimize for any cases
    def objective_funct_I_rate(self,Ir,tr,I,t,l):
        idx=np.searchsorted(t,tr)
        It=I[:,idx]
        It=It.sum(axis=0)
        I_r=Ir.sum(axis=0)
        rate=np.zeros(tr.shape[0])
        rate_t=np.zeros(tr.shape[0])
        av=np.zeros(tr.shape[0])
        av_t=np.zeros(tr.shape[0])

        rate[0]=(I_r[0]+I_r[1]-(I_r[0]-I_r[1]))/(2*I_r[0])
        rate_t[0]=(It[0]+It[1]-(It[0]-It[1]))/(2*It[0])

        for i in range(1,tr.shape[0],1):            
            rate[i]=(rate[0]-I_r[i]/I_r[i-1])/rate[0]
            rate_t[i]=(rate[0]-It[i]/It[i-1])/rate[0]
            av[i]=np.mean(rate[0:i])
            av_t[i]=np.mean(rate_t[0:i])

        return LA.norm(av_t-av)

    
    """
    # ------------------- #
    #    Miscellaneous    #
    # ------------------- #
    """

    #The tranistion model defines how to move from current to new parameters
    def transition_model(self,beta,sigma,gamma,mu,r0):
        #print("Entering transition model")
        rb=r0*(max(self.beta_r)-min(self.beta_r))
        rg=r0*(max(self.gamma_r)-min(self.gamma_r))
        rs=r0*(max(self.sigma_r)-min(self.sigma_r))
        rm=r0*(max(self.mu_r)-min(self.mu_r))
        b_p = np.random.normal(beta,rb)
        s_p = np.random.normal(sigma,rs)
        g_p = np.random.normal(gamma,rg)
        m_p = np.random.normal(mu,rm)

        while b_p > max(self.beta_r) or b_p < min(self.beta_r):
            b_p = np.random.normal(beta,rb)
        while s_p > max(self.sigma_r) or s_p < min(self.sigma_r):
            s_p = np.random.normal(sigma,rs)
        while g_p > max(self.gamma_r) or g_p < min(self.beta_r):
            g_p = np.random.normal(gamma,rg)
        while m_p > max(self.mu_r) or m_p < min(self.mu_r):
            m_p = np.random.normal(mu,rm)
        return [b_p,s_p,g_p,m_p]


    
    # Create initial points mesh
    def mesh(self,Npoints):
        print("Mesh")
        mesh = []
        for i in range(Npoints):
            beta_i=np.random.uniform(self.beta_r[0],self.beta_r[1])
            sigma_i=np.random.uniform(self.sigma_r[0],self.sigma_r[1])
            gamma_i=np.random.uniform(self.gamma_r[0],self.gamma_r[1])
            mu_i=np.random.uniform(self.mu_r[0],self.mu_r[1])
            mesh.append([beta_i,sigma_i,gamma_i,mu_i])
        return mesh


