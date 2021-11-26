from class_SEIR import SEIR
from  SEIRrefiner import SEIRrefiner 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import linalg as LA



#Import Data
if os.name == 'nt':
    path = "D:\\samue\\Dropbox\\AFES Datascience\\Ciencia y Vida\\Data\\"

else:
    path = "~/Covid_Net_SEIR/Data/"

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

dim=S_su1.shape
n=np.zeros((4*dim[0],dim[1]))
n[0:dim[0],:]=S_su1
n[dim[0]:2*dim[0],:]=E_su1
n[2*dim[0]:3*dim[0],:]=I_su1
n[3*dim[0]:4*dim[0],:]=R_su1

S0 = So[:,0]
E0 = Eo[:,0]
I0 = Io[:,0]
R0 = Ro[:,0]

#Init variables
tr=np.arange(Ir.shape[1])
h=0.1

# Parameter range
#b_r=[0.58,1.68]
b_r=[0.1,0.3] #0.2
s_r=[0.05,0.15]
#s_r=[0.05,0.15] #0.1
g_r=[0.05,0.15]
#g_r=[0.05,0.15] #0.1
mu_r=[1,3] 
#mu_r=[0.5,2.5] #2


# Strategy functions
def alpha(t):
    return(np.ones([34,34])-np.eye(34))

def eta(t):
    return(np.ones(34))

# Object init

def objective_funct(Ir,tr,I,t,l):
    idx=np.searchsorted(t,tr)
    #print(idx)
    return LA.norm(Ir-I[:,idx],l)

# SEIR Object test
mu = np.arange(1.85,2.15,0.00001)
beta = np.arange(0.185,0.215,0.00001)
sigma = np.arange(0.085,0.115,0.00001)
gamma = np.arange(0.085,0.115,0.00001)

muopt = 2
betaopt = 0.2
sigmaopt = 0.1
gammaopt = 0.1

error=[]
#Beta
# for i in range(len(beta)):
#     test = SEIR(P,eta,alpha,S0,E0,I0,R0,beta[i], sigmaopt, gammaopt,muopt)
#     test.integr(min(tr),max(tr),0.1,True)
#     error.append(objective_funct(Ir,tr,test.I,test.t,'fro'))

# print(error)
# outfile = "/home/sropert/betaerror4.csv"

# #Sigma
# for i in range(len(sigma)):
#     test = SEIR(P,eta,alpha,S0,E0,I0,R0,betaopt, sigma[i], gammaopt,muopt)
#     test.integr(min(tr),max(tr),0.1,True)
#     error.append(objective_funct(Ir,tr,test.I,test.t,'fro'))

# print(error)
# outfile = "/home/sropert/sigmaerror4.csv"

# #Gamma
# for i in range(len(gamma)):
#     test = SEIR(P,eta,alpha,S0,E0,I0,R0,betaopt, sigmaopt, gamma[i],muopt)
#     test.integr(min(tr),max(tr),0.1,True)
#     error.append(objective_funct(Ir,tr,test.I,test.t,'fro'))

# print(error)
# outfile = "/home/sropert/gammaerror4.csv"
#Mu
for i in range(len(mu)):
    test = SEIR(P,eta,alpha,S0,E0,I0,R0,betaopt, sigmaopt, gammaopt,mu[i])
    test.integr(min(tr),max(tr),0.1,True)
    error.append(objective_funct(Ir,tr,test.I,test.t,'fro'))

print(error)
outfile = "/home/sropert/muerror4.csv"


np.savetxt(outfile, error, delimiter=",") 

#error 1:
#mu = np.arange(1.8,2.2,0.002)

#error 2:
#mu = np.arange(1.98,2.02,0.0002)

#error 3:
#mu = np.arange(1.98,2.02,0.00002)

# error 4:
#mu = np.arange(1.95,2.05,0.0002)

# Seir Refiner tests
# Create param refiner object
# ref_test=SEIRrefiner(P,eta,alpha,S0,E0,I0,R0,min(tr),max(tr),0.1,b_r,g_r,s_r,mu_r)

# # # Test metropolis-hastings
# # print("Metropolis-hastings")
# # ref_test.refine(Ir,0.1,100,20,1)
# # print(ref_test.params)

# #Test  PSO
# print("PSO")
# ref_test.refinepso(n,swarmsize=50,maxiter=15,omega=0.5, phip=0.5, phig=0.5)
# # param=ref_test.paramsPSO
# # print(param)
# # d=1e-20
# ref_test2=SEIRrefiner(P,eta,alpha,S0,E0,I0,R0,min(tr),max(tr),0.1,[0.197,0.198],[0.099,0.11],[0.099,0.11],mu_r)
# ref_test2.refine(Ir,0.1,20,50,1)
# # Run integr


# # SEIR Object test
# test = SEIR(P,eta,alpha,S0,E0,I0,R0,0.21018655, 0.10047113, 0.09329693, 1.79657948)
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


# #run mesh

# #run mesh



# # Report results
