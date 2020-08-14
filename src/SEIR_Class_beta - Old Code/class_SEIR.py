#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 23:28:59 2020

@author: pmaldona
"""

import numpy as np
try:
    from scikits.odes.odeint import odeint
    scikitsimport = True
except:
    scikitsimport = False


class SEIR :

        #constructor of SEIR class elements, it's initialized when a parameter
        #miminization is performed to adjust the best setting of the actual infected

        def __init__(self,P,eta,alpha,S0,E0,I0,R0,beta,sigma,gamma,mu):
            self.scikitsimport = scikitsimport
            #init response values
            self.S=S0
            self.E=E0
            self.I=I0
            self.R=R0
            self.N=(S0+E0+I0+R0).astype("float64")
            
            #init of strategy used for the dynamics
            self.strat_prop(P,alpha,eta)
            print(self.G(1))
            #saved stragegy functions
            self.alpha=alpha
            self.eta=eta

            #saved parameter ranges for later optimization
            # init values for the coeficients
            self.beta=beta
            self.gamma=gamma
            self.sigma=sigma
            self.mu=mu
            #diferential function definitions
            self.dSdt = lambda t,S,I: -self.beta*np.diag(np.reciprocal(self.N)*S).dot(self.G(t)).dot(I);
            self.dEdt = lambda t,S,I,E: self.beta*np.diag(np.reciprocal(self.N)*S).dot(self.G(t)).dot(I) - self.sigma*E;
            self.dIdt = lambda t,E,I: self.sigma*E - self.gamma*I;
            self.dRdt = lambda t,I: self.gamma*I;

        def strat_prop(self,P,alpha,eta):

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


                    ydot[0:dim]=self.dSdt(t,y[0:dim],y[2*dim:3*dim])
                    ydot[dim:2*dim]=self.dEdt(t,y[0:dim],y[2*dim:3*dim],y[dim:2*dim])
                    ydot[2*dim:3*dim]=self.dIdt(t,y[dim:2*dim],y[2*dim:3*dim])
                    ydot[3*dim:4*dim]=self.dRdt(t,y[2*dim:3*dim]) 

                
                initcond = np.concatenate((S0, E0, I0, R0))
                # initcond = initcond.reshape(4*dim)
                soln = odeint(model_SEIR_graph, self.t, initcond)
                
                self.t=soln[1][0]   
                soln=np.transpose(np.array(soln[1][1]))
                
                self.S = soln[0:dim,:]   
                self.E = soln[dim:2*dim,:]
                self.I = soln[2*dim:3*dim,:]
                self.R = soln[3*dim:4*dim,:]

            else:
                print("Scikit couldn't be imported. Using RK4 instead")                
                return self.integr_RK4(t0,T,h,E0init)

