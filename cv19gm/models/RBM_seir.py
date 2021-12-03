
import toml
import os
import sys
#import subprocess
dirname = os.path.dirname(__file__)
import seir
#sys.path.append(dirname+"/../../KaXim-Tutorial/script")
#import plotting_kappa as kappa
sys.path.insert(1,dirname+"/../../lib")
import KaXim
import numpy as np

class RBM_SEIR(seir.SEIR):

	def kappa_sim(self,verbose = False):
		#model = toml.load(toml_file)
		model = self.cfg
		params = dict(model['parameters']['static'],**model['parameters']['dynamic'])
		params.update(model['initialconditions'])
		time = params['t_end'] - params['t_init']


		ka_params = []
		ka_params.append(params['alpha'])  #alfa
		ka_params.append(params['beta'])  #beta
		ka_params.append(params['population'])  #N
		I = params['I_det']
		if(not I):    #0 or False
			I = params['I']*params['pI_det']
		E = params['E_d']
		if(E == False):
			E = params['E']
			if(E != False):
				E = E * params['pI_det']
			else:
				E = params['mu']*I
		ka_params.append(E)  #E
		ka_params.append(I)  #I
		ka_params.append(params['R'])  #R
		ka_params.append(params['tE_I'])  #t_EI
		ka_params.append(params['tI_R'])  #t_IR
		ka_params.append(params['rR_S'])  #r_RS

		opt = ""
		for param in ka_params:
			opt = opt + " " + str(param)

		out_folder = dirname+"/runs"
		model = os.path.relpath(dirname+"/kappa/SEIR-base.xka")
		if(self.model['name'] == "SEIR-byage"):
			model = os.path.relpath(dirname+"/kappa/SEIR-by_age.xka")
		elif(self.model['name'] == "SEIR-byage-vac1"):
			model = os.path.relpath(dirname+"/kappa/SEIR-by_age-vac1.xka")
		elif(self.model['name'] == "SEIR-byage-vac2"):
			model = os.path.relpath(dirname+"/kappa/SEIR-by_age-vac2.xka")

		run_params = {
			"-i"	:	model,
			"-r"	:	"10",
			"-t"	:	str(time),
			"-p"	:	str(time),
			"--verbose"	:	"0",
			"--params"	:	" ".join([str(x) for x in ka_params])
		}
		
		if(verbose):
			print(run_params)


		#cmd_args = model+" -t "+str(time)+" -p "+str(time)+" -d "+out_folder+" --params "+opt
		#prog = subprocess.run(["KaXim",model,"-t",str(time),"-p",str(time),"-d",out_folder,"--params"]+[str(x) for x in ka_params],stdout=subprocess.PIPE,text=True)
		#if(not silence):
		#	print(prog.stdout)

		self.result = KaXim.run(run_params)
		
		sol = self.result.getAvgTrajectory().asDataFrame()
		
		self.S=sol['Susceptible'].values
		self.E=sol['Exposed'].values
		self.E_d=sol['Daily Exposed'].values
		self.I=sol['Infected'].values
		self.I_d=sol['Daily Infected'].values
		self.R=sol['Removed'].values
		self.R_d=sol['Daily Removed'].values
		#self.Flux=sol[].values

		self.E_ac = np.cumsum(self.E_d)
		#self.I_ac = np.cumsum(self.I_d) + self.I_ac
		self.R_ac = np.cumsum(self.R_d)

		#self.I_det = self.I*self.pI_det
		#self.I_d_det = self.I_d*self.pI_det
		#self.I_ac_det = self.I_ac*self.pI_det
		
		return self.result
		
	
	def solve(self,t0=0,T=None,h=0.01):
		return self.kappa_sim(false)


