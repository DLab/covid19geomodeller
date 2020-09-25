#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 23:28:59 2020

@author: pmaldona, abernardin
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
try:
    from scikits.odes.odeint import odeint
    #from scipy.integrate import odeint
    scikitsimport = True
except:
    scikitsimport = False


class SEIR:
    """SEIR with spatial nodes distribution
    
    Parameters
    ----------
    P  : array
        Mobility matrix. Shows people moving from sector i to j
    S0 : array
        Initial Susceptible array for all nodes
    E0 : array
        Initial Exposed array for all nodes
    I0 : array
        Initial Infected  array for all nodes
    R0:array
        Initial Removed array for all nodes
    beta: double
        Infection rate
    sigma: double
        Exposed rate
    gamma: double
        Removed rate

    """
    #constructor of SEIR class elements, it's initialized when a parameter
    #miminization is performed to adjust the best setting of the actual infected
    def __init__(self, P, alpha, S0, E0, I0, R0, beta, sigma, gamma, mu, epsilon, kS, kI):  # SEIR*
        self.scikitsimport = scikitsimport
        #init response values
        self.t = [0]
        self.S=S0
        self.E=E0
        self.I=I0
        self.R=R0
        self.N=(S0+E0+I0+R0).astype("float64")
        self.P=P
            
        #init of strategy used for the dynamics
        self.alpha=alpha

        #saved parameter ranges for later optimization
        # init values for the coeficients
        self.beta=beta   #infection rate
        self.gamma=gamma #removed time
        self.sigma=sigma #incubation time
        self.mu=mu       #to set initial proportion of exposed due to infected
        self.epsilon=epsilon #SEIR*
        self.kS = kS
        self.kI = kI

        ### each line should drop a 34 array ###
        self.dSdt = lambda t,S,I: -self.beta*kS*np.diag(np.reciprocal(self.N_hat(P, self.N, alpha, t))*S).dot(self.P_hat(P, alpha, t)).dot(self.I_hat(P, I, epsilon, kS, kI, alpha, t))
        self.dEdt = lambda t,S,I,E: self.beta*kS*np.diag(np.reciprocal(self.N_hat(P, self.N, alpha, t))*S).dot(self.P_hat(P, alpha, t)).dot(self.I_hat(P, I, epsilon, kS, kI, alpha, t)) - self.sigma*E;
        self.dIdt = lambda t,E,I: self.sigma*E - self.gamma*I;
        self.dRdt = lambda t,I: self.gamma*I;

    def P_hat(self, P, alpha, t): #P_hat=P(t), this implement a strategy, it modifies mobility matrix P
        return(alpha(t)*P)

    def I_hat(self, P, I, epsilon, kS, kI, alpha, t):
        return (epsilon*kS+(1-epsilon)*kI) * self.P_hat(P, alpha, t).T.dot(I)

    def N_hat(self, P, N, alpha, t):
        return(self.P_hat(P, alpha, t).T.dot(N)) # dim(N)=34
    
    def strat_prop(self,P,alpha,eta):  # Strategy propagation

        #partial definition of strategy, must be improved
        def G(t):
            return(np.diag(eta(t))+alpha(t)*P)

        self.G=G

    def integr_RK4(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.


        if(len(self.S.shape)==1):
            #pass if object is initalized
            S0=self.S
            if(E0init):
                E0=self.mu*self.I
            else:
                E0=self.E
            I0=self.I
            R0=self.R
            self.t=np.arange(t0,T+h,h)

        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition
            S0=self.S[:,idx]
            E0=self.E[:,idx]
            I0=self.I[:,idx]
            R0=self.R[:,idx]

            #set time grid
            self.t=np.arange(self.t[idx],T+h,h)


        else:
            return

        dim=self.S.shape[0]
        S=np.zeros((dim,len(self.t)))
        E=np.zeros((dim,len(self.t)))
        I=np.zeros((dim,len(self.t)))
        R=np.zeros((dim,len(self.t)))
        S[:,0]=S0
        E[:,0]=E0
        I[:,0]=I0
        R[:,0]=R0

        for j in range(len(self.t)-1):
            k1A = h*self.dSdt(self.t[j], S[:,j], I[:,j])
            k1B = h*self.dEdt(self.t[j], S[:,j], I[:,j], E[:,j])
            k1C = h*self.dIdt(self.t[j], E[:,j], I[:,j])
            k1D = h*self.dRdt(self.t[j], I[:,j])
            k2A = h*self.dSdt(self.t[j] + h/2, S[:,j] + 0.5*k1A, I[:,j] + 0.5*k1C)
            k2B = h*self.dEdt(self.t[j] + h/2, S[:,j] + 0.5*k1A, I[:,j] + 0.5*k1C, E[:,j] + 0.5*k1B)
            k2C = h*self.dIdt(self.t[j] + h/2, E[:,j] + 0.5*k1B, I[:,j] + 0.5*k1C)
            k2D = h*self.dRdt(self.t[j] + h/2, I[:,j] + 0.5*k1C)
            k3A = h*self.dSdt(self.t[j] + h/2, S[:,j] + 0.5*k2A, I[:,j] + 0.5*k2C)
            k3B = h*self.dEdt(self.t[j] + h/2, S[:,j] + 0.5*k2A, I[:,j] + 0.5*k2C, E[:,j] + 0.5*k2B)
            k3C = h*self.dIdt(self.t[j] + h/2, E[:,j] + 0.5*k2B, I[:,j] + 0.5*k2C)
            k3D = h*self.dRdt(self.t[j] + h/2, I[:,j] + 0.5*k2C)
            k4A = h*self.dSdt(self.t[j] + h, S[:,j] + k3A, I[:,j] + k3C)
            k4B = h*self.dEdt(self.t[j] + h, S[:,j] + k3A, I[:,j] + k3C, E[:,j] + k3B)
            k4C = h*self.dIdt(self.t[j] + h, E[:,j] + k3B, I[:,j] + k3C)
            k4D = h*self.dRdt(self.t[j] + h, I[:,j] + k3C)
            S[:,j+1]=S[:,j] +1/6*(k1A + 2*k2A + 2*k3A + k4A)
            E[:,j+1]=E[:,j] +1/6*(k1B + 2*k2B + 2*k3B + k4B)
            I[:,j+1]=I[:,j] +1/6*(k1C + 2*k2C + 2*k3C + k4C)
            R[:,j+1]=R[:,j] +1/6*(k1D + 2*k2D + 2*k3D + k4D)

        self.S=S
        self.E=E
        self.I=I
        self.R=R
        return
        

    def integr(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.

        if scikitsimport:
            if(len(self.S.shape)==1):
                #pass if object is initalized
                S0=self.S
                if(E0init):
                    E0=self.mu*self.I
                else:
                    E0=self.E
                I0=self.I
                R0=self.R
                self.t=np.arange(t0,T+h,h)

            elif((min(self.t)<=t0) & (t0<=max(self.t))):
                #Condition over exiting time in already initialized object

                #Search fot initial time
                idx=np.searchsorted(self.t,t0)

                #set initial condition
                S0=self.S[:,idx]
                E0=self.E[:,idx]
                I0=self.I[:,idx]
                R0=self.R[:,idx]

                #set time grid
                self.t=np.arange(self.t[idx],T+h,h)


            else:
                return()
                
            dim=self.S.shape[0]
            
            def model_SEIR_graph(t,y,ydot):
                ydot=np.zeros(4*dim)
                
                ydot[0:dim]=self.dSdt(t,y[0:dim],y[2*dim:3*dim])
                ydot[dim:2*dim]=self.dEdt(t,y[0:dim],y[2*dim:3*dim],y[dim:2*dim])
                ydot[2*dim:3*dim]=self.dIdt(t,y[dim:2*dim],y[2*dim:3*dim])
                ydot[3*dim:4*dim]=self.dRdt(t,y[2*dim:3*dim])
                
            initcond = np.concatenate((S0, E0, I0, R0))
            soln = odeint(model_SEIR_graph, self.t, initcond) #scikit

            self.t=soln[1][0]   
            soln=np.transpose(np.array(soln[1][1]))
            
            self.S = soln[0:dim,:]   
            self.E = soln[dim:2*dim,:]
            self.I = soln[2*dim:3*dim,:]
            self.R = soln[3*dim:4*dim,:]

        else:
            print("Scikit couldn't be imported. Using RK4 instead")                
            return self.integr_RK4(t0,T,h,E0init)



    def integr_sci(self,t0,T,h,E0init=False):
        #integrator function that star form t0 and finish with T with h as
        #timestep. If there aren't inital values in [t0,T] function doesn't
        #start. Or it's start if class object is initialze.
       
        if(not isinstance(self.S, np.ndarray)):
            #pass if object is initalized
            if(E0init):
                E0=self.mu*(self.I)                
            else:
                E0=self.E
                
            S0=self.S
            I0=self.I
            R0=self.R            
            
            #I_ac0=self.I_ac
            #I_d0=self.I_d

            self.t=np.arange(t0,T+h,h)
            
        elif((min(self.t)<=t0) & (t0<=max(self.t))):
            #Condition over exiting time in already initialized object

            #Search fot initial time
            idx=np.searchsorted(self.t,t0)

            #set initial condition

            S0=self.S[idx]
            E0=self.E
            
            I0=self.I[idx]
            R0=self.R[idx]                        
            #I_ac0=self.I_ac[idx]
            #I_d0=self.I_d[idx]                       

            self.t=np.arange(t0,T+h,h)

        else:
            return()
            

        
        def model_SEIR_graph(t,y):
            ydot=np.zeros(len(y))
            ydot[0]=self.dSdt(t,y[0],y[2])
            ydot[1]=self.dEdt(t,y[0],y[2],y[1])
            ydot[2]=self.dIdt(t,y[1],y[2])
            ydot[3]=self.dRdt(t,y[2])
            
            #ydot[4]=self.dI_ac(t,y[1])
            #ydot[5]=self.dI_d(t,y[1],y[5])
                         
                                          
            return(ydot)
        initcond = np.array([S0,E0,I0,R0])  

        sol = solve_ivp(model_SEIR_graph,(t0,T), initcond,method='LSODA')
        
        self.t=sol.t 

         
        self.S=sol.y[0,:]
        self.E=sol.y[1,:]
        self.I=sol.y[2,:]
        self.R=sol.y[3,:]
        #self.I_ac=sol.y[4,:]
        #self.I_d=sol.y[5,:]

        return(sol)

