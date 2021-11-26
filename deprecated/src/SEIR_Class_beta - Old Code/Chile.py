from class_SEIR import SEIR
from  SEIRrefiner import SEIRrefiner 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


#Import Data
if os.name == 'nt':
    path = "D:\\samue\\Dropbox\\AFES Datascience\\Ciencia y Vida\\Data\\"

else:
    path = "../Data/"


R_r = pd.read_excel(path+"RecupYMuerte.xlsx", header=None).to_numpy()
I_r = pd.read_excel(path+"Inf-RecupYMuerte.xlsx", header=None).to_numpy()
S0 = pd.read_excel(path+"vectorpoblacioncomuna.xlsx", header=None).to_numpy()
S0=S0.reshape(16)


tr=[-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]
tr=np.array(tr)

ind=np.where(tr[:]==7)
ind_f=np.where(tr[:]==28)

I0 = I_r[:,ind[0][0]]
R0 = R_r[:,ind[0][0]]
Ir= I_r[:,ind[0][0]:ind_f[0][0]+1]
ta=tr[ind[0][0]:ind_f[0][0]+1]
h=0.1
P=np.eye(16)
# np.where(results[:,4]==np.amin(results[:,4]))[0][0]

# Parameter range
b_r=[0.58,1.68]
#b_r=[0.1,0.3] #0.2
#s_r=[0.05,0.15]
s_r=[0.24,0.66] #0.1
#g_r=[0.05,0.15]
g_r=[0.1,0.33] #0.1
#mu_r=[1,3] 
mu_r=[0.5,2.5] #2


# Strategy functions
def alpha(t):
    return(np.zeros([16,16]))

def eta(t):
    return(np.ones(16))

# Object init


# Seir Refiner tests
# Create param refiner object
ref_test=SEIRrefiner(P,eta,alpha,S0,2*I0,I0,R0,min(ta),max(ta),h,b_r,g_r,s_r,mu_r)

# # # Test metropolis-hastings
# # print("Metropolis-hastings")
# ref_test2=SEIRrefiner(P,eta,alpha,S0,E0,I0,R0,min(tr),max(tr),0.1,[0.197,0.198],[0.099,0.11],[0.099,0.11],mu_r)
# ref_test2.refine(Ir,ta,0.1,20,50,1)
# # print(ref_test.params)

# #Test  PSO
# print("PSO")
ref_test.refinepso(Ir,ta,swarmsize=50,maxiter=15,omega=0.5, phip=0.5, phig=0.5)
param=ref_test.paramsPSO
# print(param)
# d=1e-20

# Run integr


# SEIR Object test
# test = SEIR(P,eta,alpha,S0,E0,I0,R0,0.21018655, 0.10047113, 0.09329693, 1.79657948)
# test.integr_RK4(min(tr),max(tr),0.1,True)


# mesh test

# parms=mesh(test,Ir,tr,5,5,0.025,20)
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
# plt.plot(test.t[idx],E_dif[1,:],lR_r = pd.read_excel(path+"RecupYMuerte.xlsx", header=None).to_numpy()
I_r = pd.read_excel(path+"Inf-RecupYMuerte.xlsx", header=None).to_numpy()
S0 = pd.read_excel(path+"vectorpoblacioncomuna.xlsx", header=None).to_numpy()
S0=S0[0]

tr=[-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]
tr=np.array(tr)

ind=np.where(tr[:]==7)
ind_f=np.where(tr[:]==28)

I0 = I_r[:,ind[0][0]]
R0 = R_r[:,ind[0][0]]
Ir= I_r[:,ind[0][0]:ind_f[0][0]+1]
ta=tr[ind[0][0]:ind_f[0][0]+1]
h=0.1
P=np.eye(17)
# np.where(results[:,4]==np.amin(results[:,4]))[0][0]

# Parameter range
b_r=[0.58,1.68]
#b_r=[0.1,0.3] #0.2
#s_r=[0.05,0.15]
s_r=[0.24,0.66] #0.1
#g_r=[0.05,0.15]
g_r=[0.1,0.33] #0.1
#mu_r=[1,3] 
mu_r=[0.5,2.5] #2


# Strategy functions
def alpha(t):
    return(np.zeros([17,17]))

def eta(t):
    return(np.ones(17))

# Object init


# Seir Refiner tests
# Create param refiner object
ref_test=SEIRrefiner(P,eta,alpha,S0,2*I0,I0,R0,min(ta),max(ta),h,b_r,g_r,s_r,mu_r)

# # # Test metropolis-hastings
# # print("Metropolis-hastings")
# ref_test2=SEIRrefiner(P,eta,alpha,S0,E0,I0,R0,min(tr),max(tr),0.1,[0.197,0.198],[0.099,0.11],[0.099,0.11],mu_r)
# ref_test2.refine(Ir,ta,0.1,20,50,1)
# # print(ref_test.params)

# #Test  PSO
# print("PSO")
ref_test.refinepso(Ir,ta,swarmsize=50,maxiter=15,omega=0.5, phip=0.5, phig=0.5)
param=ref_test.paramsPSO
# print(param)
# d=1e-20

# Run integr


# SEIR Object test
# test = SEIR(P,eta,alpha,S0,E0,I0,R0,0.21018655, 0.10047113, 0.09329693, 1.79657948)
# test.integr_RK4(min(tr),max(tr),0.1,True)


# mesh test

# parms=mesh(test,Ir,tr,5,5,0.025,20)
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
