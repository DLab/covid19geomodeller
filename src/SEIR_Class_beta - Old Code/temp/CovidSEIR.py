#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEIR Model param finder
Implementation of a Metropolis-Hasting model

Authors:
Pedro Maldonado, Samuel Ropert, Alejandra Barrios, Tomas Veloz
"""

import numpy as np
import pandas
from numpy import linalg as LA

class SEIRParent:
    """
    constructor of SEIR class elements, it's initialized when a parameter
    miminization is performed to adjust the best setting of the actual infected
    """
        def __init__(self,P,eta,alpha,S0,E0,I0,R0,tr):

            #init response values
            self.S=S0
            self.E=E0
            self.I=I0
            self.R=R0
            self.N=(S0+E0+I0+R0).astype("float64")

            #init of strategy used for the dynamics
            self.strat_prop(P,alpha,eta)

            #saved stragegy functions
            self.alpha=alpha
            self.eta=eta

            #definition of timegrid
            self.t0=min(tr)
            self.T=max(tr)
            self.h=h
            self.t=np.arange(self.t0,self.T+self.h,self.h)


        def strat_prop(self,P,alpha,eta):
            #partial definition of strategy, must be improved
            def G(t):
                return((np.diag(eta(t))+alpha(t))*(np.eye(P.shape[1])+P))

            self.G=G


        def integr_RK4(self,t0,T,h,beta,sigma,gamma,mu,fix=False,E0init=False):
            #integrator function that star form t0 and finish with T with h as
            #timestep. If there aren't inital values in [t0,T] function doesn't
            #start. Or it's start if class object is initialze.

            #diferential function definitions
            self.dSdt = lambda beta,t,G,N,S,I: -beta*np.diag(np.reciprocal(N)*S).dot(G(t)).dot(I);
            self.dEdt = lambda beta,sigma,t,G,N,S,I,E: beta*np.diag(np.reciprocal(N)*S).dot(G(t)).dot(I) - sigma*E;
            self.dIdt = lambda sigma,gamma,t,E,I: sigma*E - gamma*I;
            self.dRdt = lambda gamma,t,I: gamma*I;

            if(len(self.S.shape)==1):
                #pass if object is initalized
                t=self.t
                S0=self.S
                if(E0init):
                    E0=mu*self.I
                else:
                    E0=self.E
                I0=self.I
                R0=self.R

            elif((self.t0<=t0) & (t0<=self.T)):
                #Condition over exiting time in already initialized object

                #Search fot initial time
                np.where((10-h<t) & (t<10+h))
                i=np.min(np.where((10-h<t.self) & (t.self<10+h)))

                #set initial condition
                S0=self.S[:,i]
                E0=self.E[:,i]
                I0=self.I[:,i]
                R0=self.R[:,i]

                #set time grid
                t=np.arange(self.t[i],T+h,h)


            else:
                return(0)
            dim=self.S.shape[0]
            S=np.zeros((dim,len(t)))
            E=np.zeros((dim,len(t)))
            I=np.zeros((dim,len(t)))
            R=np.zeros((dim,len(t)))
            S[:,0]=S0
            E[:,0]=E0
            I[:,0]=I0
            R[:,0]=R0

            for j in range(len(t)-1):
                k1A = h*self.dSdt(beta, t[j], self.G, self.N, S[:,j], I[:,j])
                k1B = h*self.dEdt(beta, sigma, t[j], self.G, self.N, S[:,j], I[:,j], E[:,j])
                k1C = h*self.dIdt(sigma, gamma, t[j], E[:,j], I[:,j])
                k1D = h*self.dRdt(gamma, t[j], I[:,j])
                k2A = h*self.dSdt(beta, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k1A, I[:,j] + 0.5*k1C)
                k2B = h*self.dEdt(beta, sigma, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k1A, I[:,j] + 0.5*k1C, E[:,j] + 0.5*k1B)
                k2C = h*self.dIdt(sigma, gamma, t[j] + h/2, E[:,j] + 0.5*k1B, I[:,j] + 0.5*k1C)
                k2D = h*self.dRdt(gamma, t[j] + h/2, I[:,j] + 0.5*k1C)
                k3A = h*self.dSdt(beta, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k2A, I[:,j] + 0.5*k2C)
                k3B = h*self.dEdt(beta, sigma, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k2A, I[:,j] + 0.5*k2C, E[:,j] + 0.5*k2B)
                k3C = h*self.dIdt(sigma, gamma, t[j] + h/2, E[:,j] + 0.5*k2B, I[:,j] + 0.5*k2C)
                k3D = h*self.dRdt(gamma, t[j] + h/2, I[:,j] + 0.5*k2C)
                k4A = h*self.dSdt(beta, t[j] + h, self.G, self.N, S[:,j] + k3A, I[:,j] + k3C)
                k4B = h*self.dEdt(beta, sigma, t[j] + h, self.G, self.N, S[:,j] + k3A, I[:,j] + k3C, E[:,j] + k3B)
                k4C = h*self.dIdt(sigma, gamma, t[j] + h, E[:,j] + k3B, I[:,j] + k3C)
                k4D = h*self.dRdt(gamma, t[j] + h, I[:,j] + k3C)
                S[:,j+1]=S[:,j] +1/6*(k1A + 2*k2A + 2*k3A + k4A)
                E[:,j+1]=E[:,j] +1/6*(k1B + 2*k2B + 2*k3B + k4B)
                I[:,j+1]=I[:,j] +1/6*(k1C + 2*k2C + 2*k3C + k4C)
                R[:,j+1]=R[:,j] +1/6*(k1D + 2*k2D + 2*k3D + k4D)

            if(fix==True):
                self.S=S
                self.E=E
                self.I=I
                self.R=R
                self.t0=t[0]
                self.T=t[-1]
                self.t=t
                self.h=h
                self.beta=beta
                self.sigma=sigma
                self.gamma=gamma
                return(0)
            else:
                return({"S":S,"E":E,"I":I,"R":R,"t0":t[0],"T":t[-1],"h":h,
                        "t":t,"beta":beta,"sigma":sigma,"gamma":gamma})


class SEIR :
    """
    constructor of SEIR class elements, it's initialized when a parameter
    miminization is performed to adjust the best setting of the actual infected
    """
        def __init__(self,P,eta,alpha,S0,E0,I0,R0,tr,beta,gamma,sigma,mu):

            #init response values
            self.S=S0
            self.E=E0
            self.I=I0
            self.R=R0
            self.N=(S0+E0+I0+R0).astype("float64")

            #init of strategy used for the dynamics
            self.strat_prop(P,alpha,eta)

            #saved stragegy functions
            self.alpha=alpha
            self.eta=eta

            # Parameter Values
            self.beta=beta
            self.gamma=gamma
            self.sigma=sigma
            self.mu=mu

            #diferential function definitions
            self.dSdt = lambda beta,t,G,N,S,I: -beta*np.diag(np.reciprocal(N)*S).dot(G(t)).dot(I);
            self.dEdt = lambda beta,sigma,t,G,N,S,I,E: beta*np.diag(np.reciprocal(N)*S).dot(G(t)).dot(I) - sigma*E;
            self.dIdt = lambda sigma,gamma,t,E,I: sigma*E - gamma*I;
            self.dRdt = lambda gamma,t,I: gamma*I;

            #definition of timegrid
            self.t0=min(tr)
            self.T=max(tr)
            self.h=h
            self.t=np.arange(self.t0,self.T+self.h,self.h)
    """
    constructor of SEIR class elements, it's initialized when a parameter
    miminization is performed to adjust the best setting of the actual infected
    """
        def __init__(self,P,eta,alpha,S0,E0,I0,R0,tr,beta,gamma,sigma,mu):

            #init response values
            self.S=S0
            self.E=E0
            self.I=I0
            self.R=R0
            self.N=(S0+E0+I0+R0).astype("float64")

            #init of strategy used for the dynamics
            self.strat_prop(P,alpha,eta)

            #saved stragegy functions
            self.alpha=alpha
            self.eta=eta

            # Parameter Values
            self.beta=beta
            self.gamma=gamma
            self.sigma=sigma
            self.mu=mu

            #diferential function definitions
            self.dSdt = lambda beta,t,G,N,S,I: -beta*np.diag(np.reciprocal(N)*S).dot(G(t)).dot(I);
            self.dEdt = lambda beta,sigma,t,G,N,S,I,E: beta*np.diag(np.reciprocal(N)*S).dot(G(t)).dot(I) - sigma*E;
            self.dIdt = lambda sigma,gamma,t,E,I: sigma*E - gamma*I;
            self.dRdt = lambda gamma,t,I: gamma*I;

            #definition of timegrid
            self.t0=min(tr)
            self.T=max(tr)
            self.h=h
            self.t=np.arange(self.t0,self.T+self.h,self.h)


        def strat_prop(self,P,alpha,eta):
            #partial definition of strategy, must be improved
            def G(t):
                return((np.diag(eta(t))+alpha(t))*(np.eye(P.shape[1])+P))

            self.G=G


        def integr_RK4(self,t0,T,h,beta,sigma,gamma,mu,fix=False,E0init=False):
            #integrator function that star form t0 and finish with T with h as
            #timestep. If there aren't inital values in [t0,T] function doesn't
            #start. Or it's start if class object is initialze.


            if(len(self.S.shape)==1):
                #pass if object is initalized
                t=self.t
                S0=self.S
                if(E0init):
                    E0=mu*self.I
                else:
                    E0=self.E
                I0=self.I
                R0=self.R

            elif((self.t0<=t0) & (t0<=self.T)):
                #Condition over exiting time in already initialized object

                #Search fot initial time
                np.where((10-h<t) & (t<10+h))
                i=np.min(np.where((10-h<t.self) & (t.self<10+h)))

                #set initial condition
                S0=self.S[:,i]
                E0=self.E[:,i]
                I0=self.I[:,i]
                R0=self.R[:,i]

                #set time grid
                t=np.arange(self.t[i],T+h,h)


            else:
                return(0)
            dim=self.S.shape[0]
            S=np.zeros((dim,len(t)))
            E=np.zeros((dim,len(t)))
            I=np.zeros((dim,len(t)))
            R=np.zeros((dim,len(t)))
            S[:,0]=S0
            E[:,0]=E0
            I[:,0]=I0
            R[:,0]=R0

            for j in range(len(t)-1):
                k1A = h*self.dSdt(beta, t[j], self.G, self.N, S[:,j], I[:,j])
                k1B = h*self.dEdt(beta, sigma, t[j], self.G, self.N, S[:,j], I[:,j], E[:,j])
                k1C = h*self.dIdt(sigma, gamma, t[j], E[:,j], I[:,j])
                k1D = h*self.dRdt(gamma, t[j], I[:,j])
                k2A = h*self.dSdt(beta, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k1A, I[:,j] + 0.5*k1C)
                k2B = h*self.dEdt(beta, sigma, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k1A, I[:,j] + 0.5*k1C, E[:,j] + 0.5*k1B)
                k2C = h*self.dIdt(sigma, gamma, t[j] + h/2, E[:,j] + 0.5*k1B, I[:,j] + 0.5*k1C)
                k2D = h*self.dRdt(gamma, t[j] + h/2, I[:,j] + 0.5*k1C)
                k3A = h*self.dSdt(beta, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k2A, I[:,j] + 0.5*k2C)
                k3B = h*self.dEdt(beta, sigma, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k2A, I[:,j] + 0.5*k2C, E[:,j] + 0.5*k2B)
                k3C = h*self.dIdt(sigma, gamma, t[j] + h/2, E[:,j] + 0.5*k2B, I[:,j] + 0.5*k2C)
                k3D = h*self.dRdt(gamma, t[j] + h/2, I[:,j] + 0.5*k2C)
                k4A = h*self.dSdt(beta, t[j] + h, self.G, self.N, S[:,j] + k3A, I[:,j] + k3C)
                k4B = h*self.dEdt(beta, sigma, t[j] + h, self.G, self.N, S[:,j] + k3A, I[:,j] + k3C, E[:,j] + k3B)
                k4C = h*self.dIdt(sigma, gamma, t[j] + h, E[:,j] + k3B, I[:,j] + k3C)
                k4D = h*self.dRdt(gamma, t[j] + h, I[:,j] + k3C)
                S[:,j+1]=S[:,j] +1/6*(k1A + 2*k2A + 2*k3A + k4A)
                E[:,j+1]=E[:,j] +1/6*(k1B + 2*k2B + 2*k3B + k4B)
                I[:,j+1]=I[:,j] +1/6*(k1C + 2*k2C + 2*k3C + k4C)
                R[:,j+1]=R[:,j] +1/6*(k1D + 2*k2D + 2*k3D + k4D)

            if(fix==True):
                self.S=S
                self.E=E
                self.I=I
                self.R=R
                self.t0=t[0]
                self.T=t[-1]
                self.t=t
                self.h=h
                self.beta=beta
                self.sigma=sigma
                self.gamma=gamma
                return(0)
            else:
                return({"S":S,"E":E,"I":I,"R":R,"t0":t[0],"T":t[-1],"h":h,
                        "t":t,"beta":beta,"sigma":sigma,"gamma":gamma})


        def strat_prop(self,P,alpha,eta):
            #partial definition of strategy, must be improved
            def G(t):
                return((np.diag(eta(t))+alpha(t))*(np.eye(P.shape[1])+P))

            self.G=G


        def integr_RK4(self,t0,T,h,beta,sigma,gamma,mu,fix=False,E0init=False):
            #integrator function that star form t0 and finish with T with h as
            #timestep. If there aren't inital values in [t0,T] function doesn't
            #start. Or it's start if class object is initialze.


            if(len(self.S.shape)==1):
                #pass if object is initalized
                t=self.t
                S0=self.S
                if(E0init):
                    E0=mu*self.I
                else:
                    E0=self.E
                I0=self.I
                R0=self.R

            elif((self.t0<=t0) & (t0<=self.T)):
                #Condition over exiting time in already initialized object

                #Search fot initial time
                np.where((10-h<t) & (t<10+h))
                i=np.min(np.where((10-h<t.self) & (t.self<10+h)))

                #set initial condition
                S0=self.S[:,i]
                E0=self.E[:,i]
                I0=self.I[:,i]
                R0=self.R[:,i]

                #set time grid
                t=np.arange(self.t[i],T+h,h)


            else:
                return(0)
            dim=self.S.shape[0]
            S=np.zeros((dim,len(t)))
            E=np.zeros((dim,len(t)))
            I=np.zeros((dim,len(t)))
            R=np.zeros((dim,len(t)))
            S[:,0]=S0
            E[:,0]=E0
            I[:,0]=I0
            R[:,0]=R0

            for j in range(len(t)-1):
                k1A = h*self.dSdt(beta, t[j], self.G, self.N, S[:,j], I[:,j])
                k1B = h*self.dEdt(beta, sigma, t[j], self.G, self.N, S[:,j], I[:,j], E[:,j])
                k1C = h*self.dIdt(sigma, gamma, t[j], E[:,j], I[:,j])
                k1D = h*self.dRdt(gamma, t[j], I[:,j])
                k2A = h*self.dSdt(beta, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k1A, I[:,j] + 0.5*k1C)
                k2B = h*self.dEdt(beta, sigma, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k1A, I[:,j] + 0.5*k1C, E[:,j] + 0.5*k1B)
                k2C = h*self.dIdt(sigma, gamma, t[j] + h/2, E[:,j] + 0.5*k1B, I[:,j] + 0.5*k1C)
                k2D = h*self.dRdt(gamma, t[j] + h/2, I[:,j] + 0.5*k1C)
                k3A = h*self.dSdt(beta, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k2A, I[:,j] + 0.5*k2C)
                k3B = h*self.dEdt(beta, sigma, t[j] + h/2, self.G, self.N, S[:,j] + 0.5*k2A, I[:,j] + 0.5*k2C, E[:,j] + 0.5*k2B)
                k3C = h*self.dIdt(sigma, gamma, t[j] + h/2, E[:,j] + 0.5*k2B, I[:,j] + 0.5*k2C)
                k3D = h*self.dRdt(gamma, t[j] + h/2, I[:,j] + 0.5*k2C)
                k4A = h*self.dSdt(beta, t[j] + h, self.G, self.N, S[:,j] + k3A, I[:,j] + k3C)
                k4B = h*self.dEdt(beta, sigma, t[j] + h, self.G, self.N, S[:,j] + k3A, I[:,j] + k3C, E[:,j] + k3B)
                k4C = h*self.dIdt(sigma, gamma, t[j] + h, E[:,j] + k3B, I[:,j] + k3C)
                k4D = h*self.dRdt(gamma, t[j] + h, I[:,j] + k3C)
                S[:,j+1]=S[:,j] +1/6*(k1A + 2*k2A + 2*k3A + k4A)
                E[:,j+1]=E[:,j] +1/6*(k1B + 2*k2B + 2*k3B + k4B)
                I[:,j+1]=I[:,j] +1/6*(k1C + 2*k2C + 2*k3C + k4C)
                R[:,j+1]=R[:,j] +1/6*(k1D + 2*k2D + 2*k3D + k4D)

            if(fix==True):
                self.S=S
                self.E=E
                self.I=I
                self.R=R
                self.t0=t[0]
                self.T=t[-1]
                self.t=t
                self.h=h
                self.beta=beta
                self.sigma=sigma
                self.gamma=gamma
                return(0)
            else:
                return({"S":S,"E":E,"I":I,"R":R,"t0":t[0],"T":t[-1],"h":h,
                        "t":t,"beta":beta,"sigma":sigma,"gamma":gamma})



def mesh(Npoints):
    mesh = []
    for i in range(Npoints):
        beta_i=np.random.uniform(model.beta_r)
        sigma_i=np.random.uniform(model.sigma_r)
        gamma_i=np.random.uniform(model.gamma_r)
        mu_i=np.random.uniform(model.mu_r)
        mesh.append([beta_i,sigma_i,gamma_i,mu_i])
    return mesh

## Definir un objeto SEIR que solo se inicialice con las variables que no cambian
## Luego definir uno heredado que se defina con esas variables y los parametros
def met_hast(model,beta_i,sigma_i,gamma_i,mu_i,r0,steps,err):
    x=model.integr_RK4(model.t0,model.T,model.h,beta_i,sigma_i,gamma_i,mu_i,False,True)
    e_0=objective_funct(model.Ir,model.tr,x["I"],x["t"])
    e_o=e_0
    params = [[beta_i,sigma_i,gamma_i,mu_i,e_o]]
    for i in range(steps):
        [b_p,s_p,g_p,m_p]=transition_model(x["beta"],x["sigma"],x["gamma"],x["mu"],r0,model)
        x_new=model.integr_RK4(model.t0,model.T,model.h,b_p,s_p,g_p,m_p,False,True)
        e_n=objective_funct(model.Ir,model.tr,x_new["I"],x_new["t"])
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
    rb=r0*(max(model.beta_r)-min(model.beta_r))
    rg=r0*(max(model.gamma_r)-min(model.gamma_r))
    rs=r0*(max(model.sigma_r)-min(model.sigma_r))
    rm=r0*(max(model.mu_r)-min(model.mu_r))
    b_p = np.random.normal(beta,rb)
    s_p = np.random.normal(sigma,rs)
    r_p = np.random.normal(gamma,rg)
    m_p = np.random.normal(mu_i,rm)

    while b_p > max(model.beta_r) or b_p < min(model.beta_r):
        b_p = np.random.normal(beta,rb)
    while s_p > max(model.sigma_r) or s_p < min(model.sigma_r):
        s_p = np.random.normal(sigma,rs)
    while g_p > max(model.gamma_r) or g_p < min(model.beta_r):
        g_p = np.random.normal(gamma,rg)
    while m_p > max(model.mu_r) or m_p < min(model.mu_r):
        m_p = np.random.normal(beta,rb)

    return [b_p,s_p,g_p,m_p]

# Strategy functions
def alpha(t):
    return(np.ones([34,34])-np.eye(34))

def eta(t):
    return(np.ones(34))


if __name__ == "__main__":
    #Import Data
    So = pd.read_excel("poblacion_Inicial_S_stgo.xlsx", header=None).to_numpy()
    S0 = So[:,0]
    Eo = pd.read_excel("poblacion_Inicial_E_stgo.xlsx", header=None).to_numpy()
    E0 = Eo[:,0]
    Io = pd.read_excel("poblacion_Inicial_I_stgo.xlsx", header=None).to_numpy()
    I0 = Io[:,0]
    Ro = pd.read_excel("poblacion_Inicial_R_stgo.xlsx", header=None).to_numpy()
    R0 = Ro[:,0]
    P = pd.read_excel("connectivity_stgo2.xlsx", header=None).to_numpy()
    Ir=pd.read_excel("Simulacion-400dias-I.xlsx", header=None).to_numpy()


    #Init variables
    tr=np.arange(Ir.shape[1])
    h=0.1
    r0 = 0.1
    steps = 20
    err = 20

    # Parameter range
    b_r=[0.01,0.5] #0.1
    s_r=[0.01,0.5] #0.1
    g_r=[0.01,0.5] #0.1
    mu_r=[0.5,4] #2

    SEIRobject = SEIRParent(P,eta(),alpha(),S0,E0,I0,R0,tr)

    #Create Mesh
    Npoints = 20 # Number of initial points in params hyperspace
    mesh = mesh(Npoints)
    results = []

    for i in range(Npoints):
        results.append(met_hast(SEIRobject,mesh[0],mesh[1],mesh[2],mesh[3],r0,steps,err)[-1])

    # Find the best params:
    min(results[:,3])
    results.index(min(results[:,3,]))


    # init model
    #params = [[beta_i,sigma_i,gamma_i,mu_i,e_o]]
    #model = SEIR(self,P,eta,alpha,S0,E0,I0,R0,Ir,tr,h,beta_r,gamma_r,sigma_r,mu_r)
    # Le sacaria vairables de inicializacion al modelo SeIR
    # no se usan beta, gamma, sigma, mu,y h
    # I_r se usa en todos los modelos, lo abstraeria un nivel
    #   def __init__(self,P,eta,alpha,S0,E0,I0,R0,Ir,tr,h,beta_r,gamma_r,sigma_r,mu_r):
#           n = 5 # Number of initial parameters
    # create a mesh of initial parameters


    # execute met_hast for each parameter

    # Output comparison
