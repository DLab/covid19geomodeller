#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#from logging import raiseExceptions
from distutils.log import error
import matplotlib.pyplot as plt
import numpy as np
import time

# optimization
import pygmo as pg
from numpy import linalg as LA

# cv19gm libraries
from cv19gm.data.cv19data import ImportData as importdata
from cv19gm.cv19sim import CV19SIM
import cv19gm.utils.cv19functions as cv19functions
from cv19gm.utils.cv19errors import *

"""To Do
* Implementar método Alejo:
    * Ver si se puede optimizar, demora como media hora =S 

* Implementar mi método
* Agregarlo a la función general cv19paramfit


# To Do:
[x] 1. Arreglar lo del cálculo de error para datos discontinuos
[x] 2. Testear error ScaledRMSE vs RMSE
[x] 3. Crear función piecewise con sólo los días de intermedios [t0,t1,t2,t3] en vez de [[t0,t1],[t2,t3],[t4,t5]]
[x] 4. Comprobar que el opti_mu esté bien implementado
[x] 5. Intentar no tener que simular nuevamente con el beta obtenido, ver si se puede obtener la última simulación: can't
[ ] 6. Calcular error global al avanzar con las iteraciones
[ ] 7. Implementar número mínimo de días entre cambios de tendencias
[ ] 8. Implementar iteración beta_fit. Condiciones de salida:
       * Error < erro_tol
       * Numero de iteraciones
       * Disminusión relativa de error entre iteraciones
[ ] 9. Revisar consistencia al establecer las condiciones iniciales de las simulaciones
[ ] 10. Agregar la posibilidad de meter objetos de datos a los fiteos
[ ] 11. Ordenar, limpiar y reducir el código
[ ] 12. Exponer variables de optimización como población y número de generaciones

"""


class CV19PARAMFIT:
    def __init__(
        self,
        cfg,
        method="interval_fit",
        I_d_data=None,
        t_data=None,
        dates_data=None,
        tstate=None,
        initdate=None,
        cv19gmdata=None,
        verbose=False,
        nintervals=3,
        gen=200,
        vartofit="I_d",
        **kwargs
    ):
        """CV19PARAMFIT
        Data fitting library for cv19gm platform. Manages the data and initializes the optimization method.

        The data to fit can be given in three different ways:
            * Data array using I_d_data and t_data
            * Target state id and initial date for obtaining the data from our servers, usubg tstate and initdate
            * cv19gmdata object with the data already loaded

        Args:
            cfg (str,dict): Simulation configuration file
            method (str, optional): Fitting method (add list of available methods). Defaults to 'interval_fit'.
            I_d_data (_type_, optional): Data new cases to fit. Defaults to None.
            t_data (list or nparray, optional): Data time array. Defaults to None.
            dates_data (_type_, optional): Data dates array. Defaults to None.
            tstate (_type_, optional): State code to obtain data to fit. Defaults to None.
            initdate (_type_, optional): Initial date to fit from data obtained automatically. Defaults to None.
            cv19gmdata (_type_, optional): cv19gmdata object with data to fit. Defaults to None.
            verbose (bool, optional): Verbose output. Defaults to False.
            nintervals (int, optional): Number of intervals to fit. Applies to interval_fit method. Defaults to 3.
            gen (int, optional): Number of generation for fitting method. Defaults to 200.
            vartofit (str, optional): Variable to fit on simulation. Does nothing for the moment, but it may be used in the future. Defaults to 'I_d'.
        """
        self.cfg = cfg
        self.gen = gen
        self.kwargs = kwargs

        # Data management: data source and data load

        # Data array
        if not type(I_d_data) == type(None):
            self.I_d_data = I_d_data
            self.t_data = t_data
            self.dates_data = dates_data
            self.data = None

        # State code + inital date: downloads data from remote server
        elif not type(tstate) == type(None):
            self.data = importdata(tstate=tstate, initdate=initdate)
            # Import data for seir - change this for adding more models
            print("Importing data from remote server")
            self.data.import_data_seir()
            self.I_d_data, self.t_data, self.dates_data = (
                self.data.I_d_r,
                self.data.I_d_r_tr,
                self.data.I_d_r_dates,
            )

        # cv19gmdata object
        elif not type(cv19gmdata) == type(None):
            self.data = cv19gmdata
            self.I_d_data, self.t_data, self.dates_data = (
                self.data.I_d_r,
                self.data.I_d_r_tr,
                self.data.I_d_r_dates,
            )

        # Missing data to fit
        else:
            print("Missing data to fit")
            return
            # raiseExceptions

        # Optimization methods

        if method == "interval_fit":
            self.nintervals = nintervals  # number of intervals

            # Set bounds
            alpha_low_bound = np.repeat(0.01, self.nintervals)  # lower bound for alpha
            days_interval_low_bound = np.repeat(
                0, self.nintervals
            )  # lower bound for days interval

            alpha_up_bound = np.repeat(1.0, self.nintervals)  # upper bound alpha
            days_interval_up_bound = np.repeat(
                90, self.nintervals
            )  # upper bound for days interval, 90 is the max day between intervals.

            low_bound = np.append(alpha_low_bound, days_interval_low_bound)
            up_bound = np.append(alpha_up_bound, days_interval_up_bound)

            bounds = [low_bound, up_bound]
            t_sim = 500  # simulation time in days

            # until = 150                                   # set limit to optimize with real data

            # Set individual
            # ver si funciona bien la inicialización de I_d para casos sin inputdata
            self.opti = INTERVAL_FIT(
                I_d_data=self.I_d_data,
                t_data=self.t_data,
                bounds=bounds,
                cfg=cfg,
                nintervals=self.nintervals,
                inputdata=self.data,
                I_d=self.I_d_data[0],
            )

            # Set Particle Swarm Optimization parameters, gen=200, population=40 -> time: 20 min
            self.algo = pg.algorithm(pg.pso(gen=self.gen))
            self.algo.set_verbosity(10)

        elif method == "other_method":
            print("Performing other method =)")

    def evolve(self, npop=40):
        self.pop = pg.population(self.opti, npop)
        self.pop = self.algo.evolve(self.pop)
        self.opti.fitness(self.pop.champion_x)
        print(self.pop.champion_f)

    def simulate(self):
        if not type(self.dates_data) == type(None):
            self.kwargs["initdate"] = self.dates_data[0]
        self.sim = CV19SIM(
            self.cfg,
            beta=self.opti.beta_funct,
            inputdata=self.data,
            I_d=self.I_d_data[0],
            **self.kwargs
        )
        self.sim.solve()

    def plot(self, xaxis="days"):
        """Plot simulation results with optimized parameters"""
        interval_array = [
            self.opti.beta_funct["days"][i][1]
            for i in range(len(self.opti.beta_funct["days"]))
        ]
        if xaxis == "days":
            plt.plot(self.sim.t, self.sim.I_d, label="Simulation")
            plt.scatter(self.t_data, self.I_d_data, label="Data")
            plt.xlabel("Days")
            for i in np.array(interval_array):
                plt.axvline(i, linestyle="dashed", color="grey")
        elif xaxis == "dates":
            plt.plot(self.sim.dates, self.sim.I_d, label="Simulation")
            plt.scatter(self.dates_data, self.I_d_data, label="Data")
            plt.xlabel("Dates")
            for i in np.array(interval_array):
                plt.axvline(
                    self.sim.dates[int(i)], linestyle="dashed", color="grey"
                )

        plt.ylabel("New cases")
        plt.title("Optimization results")

        # plot vertical lines, where days intervals start

    def plot_history(self):
        """Plots optimization history
        # log of the optimization, for PSO. log = [Gen, Fevals, gbest, Mean Vel., Mean lbest, Avg. Dist]
        # Gen: generation, Fevals:fitness evaluations, gbest: global best,
        # Mean Vel.: mean velocity, Mean lbest: mean fitness, Avg. Dist: average distance
        """
        uda = self.algo.extract(pg.pso)  # uda = user-defined algorithm
        log = np.array(uda.get_log())
        # plot global best vs generations
        plt.plot(log[:, 0], log[:, 2])
        plt.xlabel("Generation", size=15)
        plt.ylabel("Fitness", size=15)
        plt.show()
        return


class BETA_FIT:
    """
    Optimizador de beta con Inverse Time Weighted Error


    1. Se deben encontrar las condiciones iniciales para todas las variables
    2. Elegir qué variable se utilizará en el fitness. Voy a partir con I_d


    """

    def __init__(self, I_d, t, cfg, bounds, error="RMSE", alpha=1, **kwargs):
        self.I_d = I_d  # New daily infected
        self.t = t  # real time
        self.bounds = bounds  # Bounds for the variables to optimize
        self.cfg = cfg
        self.alpha = alpha

        self.kwargs = kwargs

        if error == "ITWE":
            self.error = ITWE
        elif error == "RMSE":
            self.error = RMSE

        self.ITWE_error = 0
        self.RMSE_error = 0
        self.p2p_error = []

        #
        # self.mu = mu

    def fitness(self, x):
        sim = CV19SIM(self.cfg, beta=x, **self.kwargs)
        sim.solve()
        self.ITWE_error = ITWE(sim.I_d, self.I_d, self.t, alpha=self.alpha)
        self.RMSE_error = RMSE(sim.I_d, self.I_d, self.t)
        self.p2p_error = P2PE(sim.I_d, self.I_d, self.t)
        return [self.ITWE_error]

    def get_bounds(self):  # mandatory function of Pygmo2
        return self.bounds

    def set_bounds(self, bounds):
        self.bounds = bounds
        return self.bounds

    def gradient(self, x):
        return pg.estimate_gradient_h(lambda x: self.fitness(x), x)


class BETAMU_FIT:
    """
    Optimizador de beta con Inverse Time Weighted Error

    1. Se deben encontrar las condiciones iniciales para todas las variables
    2. Elegir qué variable se utilizará en el fitness. Voy a partir con I_d
    """
    def __init__(self, cfg, bounds, I_d_data, t_data = None,  error="ITWE", rho=1, **kwargs):
       
        self.I_d_data = I_d_data  # New cases 
        self.t_data = list(range(len(I_d_data))) if type(t_data) == type(None) else t_data
        self.bounds = bounds  # Bounds for the variables to optimize
        self.cfg = cfg
        self.rho = rho # Decide if keeping this
        self.kwargs = kwargs
        if error == "ITWE":
            self.error = ITWE
        elif error == "RMSE":
            self.error = RMSE

    def fitness(self, x):  # mandatory function of Pygmo2
        # We construct the days intervals with the values that are being
        # optimized.
        sim = CV19SIM(self.cfg, beta=x[0], mu=x[1], I_d = self.I_d_data[0],**self.kwargs)
        sim.solve()
        err = self.error(sim.I_d, self.I_d_data, t_data=self.t_data)
        return [err]

    def get_bounds(self):  # mandatory function of Pygmo2
        return self.bounds

    def set_bounds(self, bounds):
        self.bounds = bounds
        return self.bounds

    def gradient(self, x):
        return pg.estimate_gradient_h(lambda x: self.fitness(x), x)

class BETA_PIECEWISE_FIT:
    """_summary_
    # beta_days includes the last transition date, the objetive of this method is to find the best value for this transition.
    """
    def __init__(self,cfg, bounds, I_d_data, t_data = None,  beta_values = [], beta_days = [-np.inf], error="RMSE", **kwargs):
        
        self.I_d_data = I_d_data  # Data new cases 
        if type(t_data) == type(None):
            t_data = list(range(len(I_d_data)))
        self.t_data = t_data  # Data time
        self.bounds = bounds  # Bounds for the variables to optimize [[beta_min],[beta_max]]
        self.cfg = cfg # Configuration file for the rest of the model's parameters
        self.beta_values = beta_values # if none we start optimizing from the beginning
        self.beta_days = beta_days 
        self.kwargs = kwargs

        if error == "ITWE":
            self.error = ITWE
        elif error == "RMSE":
            self.error = RMSE
        elif error == "RRMSE":
            self.error = RRMSE
        else:
            raise Exception('Wrong error function')

    def fitness(self, x):
        beta = cv19functions.piecewise(values=self.beta_values+[x],limits = self.beta_days)
        sim = CV19SIM(self.cfg, beta=beta, I_d = self.I_d_data[0], **self.kwargs)
        sim.solve()
        self.RMSE_error = RMSE(sim.I_d, self.I_d_data, self.t_data)
        return [self.RMSE_error]

    def get_bounds(self):  # mandatory function of Pygmo2
        return self.bounds

    def set_bounds(self, bounds):
        self.bounds = bounds
        return self.bounds

    def gradient(self, x):
        return pg.estimate_gradient_h(lambda x: self.fitness(x), x)

class INTERVAL_FIT:
    """Data Fit
    Developed by A. Bernardin, modified by S. Ropert
    """
    def __init__(
        self, I_d_data, t_data, bounds, cfg, error="RMSE", nintervals=10, **kwargs
    ):

        self.I_d_data = I_d_data  # real infected
        self.t_data = t_data  # real time
        self.bounds = bounds
        self.alpha_interval_array = 0  # interval_array
        self.beta_funct = ""
        self.cfg = cfg
        self.kwargs = kwargs
        self.nintervals = nintervals

        # if error == 'ITWE':
        #     self.error =  ITWE
        # elif error == 'RMSE':
        #     self.error =  RMSE

    def fitness(self, x):  # mandatory function of Pygmo2
        """
        We search a vector that includes values for alpha and for the days intervals.

        x: Optimization vector, by default from Pygmo2. Vector of lenght 20, where x[:10] are the alpha values and x[10:] are
           the days values for the intervals.
        """
        # We construct the days intervals with the values that are being
        # optimized.
        alpha_interval_array = []
        first_element = 0
        for i in x[self.nintervals :]:
            second_element = first_element + i
            alpha_interval_array.append([first_element, second_element])
            first_element = second_element

        self.alpha_interval_array = alpha_interval_array
        beta_funct = {
            "function": "events",
            "values": [],
            "days": [],
        }  # dictionary that includes betas and values for days intervals
        beta_funct["values"] = list(x[: self.nintervals])
        beta_funct["days"] = alpha_interval_array

        self.beta_funct = beta_funct

        # set simulation
        sims = CV19SIM(
            self.cfg, beta=beta_funct, t_end=self.t_data[-1] + 1, **self.kwargs
        )  # x[i] takes values between limits define by "bounds"
        sims.solve()  # run simulation
        idx = np.searchsorted(
            sims.t, self.t_data
        )  # to adjust data to simulation (data are each three days)
        res = LA.norm(
            self.I_d_data - sims.I_d[idx]
        )  # Norm between data and simulation. This is the fitness.
        return [res]

    def get_bounds(self):  # mandatory function of Pygmo2
        return self.bounds

    def set_bounds(self, bounds):
        self.bounds = bounds
        return self.bounds


class SEQUENTIAL_FIT:
    """
    Sequential fit optimization method
    1. We aply betamu_fit for fiting the initial parameters using a truncated part of the data
    2. We calculate a global error function and check if we are under the error tolerance. If not we continue with the following points
    3. Using tendency_change we look for the tendency change day.
    4. We optimize with beta_fit from the day found in (2). 
    5. We iterate points 2, 3 and 4 until we reach the end of the data or we reduce the error below the tolerance

    """
    def __init__(self, cfg, I_d_data=None, t_data=None, dates_data=None, tstate=None,
        initdate=None, cv19gmdata=None, global_errortol=100, local_errortol=100, global_errfunct=RMSE,local_errfunct=LAE, 
        intervalsize=10, maxintervals=5, bounds_beta=[0,1], bounds_mu=[0,4], inputdata = None, paramtol=0.05, **kwargs):

        self.cfg = cfg
        self.I_d_data = I_d_data
        self.t_data = t_data
        self.dates_data = dates_data
        self.tstate = tstate
        self.initdate = initdate
        self.cv19gmdata = cv19gmdata

        self.cfg = cfg
        self.inputdata = inputdata
        self.kwargs = kwargs
        self.intervalsize = intervalsize
        self.maxintervals = maxintervals

        self.bounds_beta = bounds_beta
        self.bounds_mu = bounds_mu

        self.mu = 1
        self.beta_values = []
        self.beta_days = []

        self.global_errortol = global_errortol
        self.local_errortol = local_errortol 
        
        self.stop = False

        self.actualidx = intervalsize
        self.paramtol = paramtol
        
        # Data management: data source
        data_management(self, I_d_data, t_data, dates_data, tstate, initdate, cv19gmdata)

        # Choose error functions (later will add the capability of choosing)
        self.global_errfunct = cv19errorbuild(global_errfunct)
        self.local_errfunct = cv19errorbuild(local_errfunct)
        
        self.betaiter = 0
        self.global_error_hist = []

    # ------------------------- #
    # First step: betamu_fit    #
    # ------------------------- #
    def betamu_fit(self,actualidx=False, debug = True):
        if not actualidx:
            actualidx = self.intervalsize
        start = time.time()
        # Set bounds
        lb = [self.bounds_beta[0], self.bounds_mu[0]]  # lower bound
        ub = [self.bounds_beta[1], self.bounds_mu[1]]  # upper bound
        bounds = [lb, ub]

        # Build beta_mu optimization Problem
        opti_mu = BETAMU_FIT(cfg = self.cfg, bounds = bounds, I_d_data = self.I_d_data[:actualidx], t_data = self.t_data[:actualidx], error="RMSE", inputdata = self.inputdata, **self.kwargs)
        # Choose algorithm
        algo = pg.algorithm(pg.nlopt(solver="bobyqa"))
        algo.set_verbosity(10)

        # Solve
        if debug == True:
            print("Solving for beta and mu")
        pop = pg.population(opti_mu, size=20)
        pop = algo.evolve(pop)

        # Extract results
        self.beta_values.append(np.around(pop.champion_x[0],4))
        self.mu = np.around(pop.champion_x[1],4)

        # Calculate Global Error
        # Try to avoid simulating again
        self.sim = CV19SIM(self.cfg,I_d = self.I_d_data[0],beta=self.beta_values[0],mu=self.mu,inputdata = self.inputdata, **self.kwargs)
        self.sim.solve()
        self.global_error = self.global_errfunct(self.sim.I_d, self.I_d_data, t_data=self.t_data)
        self.global_error_hist.append(self.global_error)

        end = time.time()
        if debug == True:
            print('Initial values')
            print('Elapsed time:'+str(np.around(end - start,2))+' s')
            print('beta: '+str(self.beta_values[0]))
            print('mu '+str(self.mu))  
            print('Global Error: '+str(self.global_error))  

    # ----------------------------------------------------- #
    # Second step: recursing beta fit for short intervals    #
    # ----------------------------------------------------- #
    def beta_fit(self,actualidx=-1,debug=False):
        starttime = time.time()
        # Find right beta for the divergence point
        bounds = [[self.bounds_beta[0]],[self.bounds_beta[1]]]
        opti = BETA_PIECEWISE_FIT(cfg = self.cfg, bounds = bounds, I_d_data = self.I_d_data[:actualidx], 
            t_data = self.t_data[:actualidx],  beta_values=self.beta_values,beta_days=self.beta_days, inputdata = self.inputdata, mu = self.mu,  **self.kwargs)

        # Choose algorithm
        algo = pg.algorithm(pg.nlopt(solver="bobyqa"))
        algo.set_verbosity(10)

        # Find parameters
        if debug:
            print("Solving for beta")
        pop = pg.population(opti, size=20)
        pop = algo.evolve(pop)

        self.beta_values.append(np.around(pop.champion_x[0],4))

        # Calculate Global error
        self.beta = cv19functions.piecewise(values=self.beta_values,limits=self.beta_days)
        self.sim = CV19SIM(self.cfg,I_d = self.I_d_data[0],beta=self.beta,mu=self.mu,inputdata=self.inputdata,**self.kwargs)
        self.sim.solve()

        self.global_error = self.global_errfunct(self.sim.I_d, self.I_d_data, t_data=self.t_data)
        self.global_error_hist.append(self.global_error)
        
        endtime = time.time()
        
        if debug:
            print('Optimization process Finished')
            print('Elapsed time:'+str(np.around(endtime - starttime,2))+' s')
            print('beta: '+str(self.beta_values))
            print('days: '+str(self.beta_days))
            print('Global Error: '+str(self.global_error))            
        
        self.betaiter+=1

    def optimize(self,debug=False):
        # Finding initial values
        start = time.time()
        self.betamu_fit(debug=debug)
        
        # Stop condition: global error under a tolerance or reaching the end of the data
        while self.global_error > self.global_errortol and not self.stop:
            
            # 1. Find divergence point 
            idx = tendency_change(self.I_d_data,  t_data=self.t_data, sim = self.sim, errorfunction = self.local_errfunct, l_err_tol=self.local_errortol,startingday=self.actualidx) - 1
            
            # Choose the last data point to fit in this iteration. Check if we reached the end of the data
            if idx+self.intervalsize > len(self.I_d_data):
                if debug:
                    print('Reached the end of data')
                self.actualidx = len(self.I_d_data) 
                self.stop = True
            else:
                self.actualidx = idx+self.intervalsize

            # 2. Find the best fit for this interval
            self.beta_days.append(idx)
            self.beta_fit(debug=debug,actualidx=self.actualidx)
            
            # 3. Check if the solution represents a transition or it is due to noise.
            #    If the betas are too similar we find one that fits this extended range altogether
            if np.abs(self.beta_values[-1] - self.beta_values[-2])<self.paramtol:
                # Delete last idx and 2 last betas
                del self.beta_days[-1]
                del self.beta_values[-2:]
                # Find beta for the extended range
                if len(self.beta_days)>0:
                    self.beta_fit(debug=False,actualidx=self.actualidx)
                else:
                    self.betamu_fit(debug=False,actualidx=self.actualidx)

        end = time.time()

        print('Optimization process Finished')
        print('Elapsed time:'+str(np.around(end - start,2))+' s')
        print('beta: '+str(self.beta_values))
        print('days: '+str(self.beta_days))
        print('mu: '+str(self.mu))
        print('Global Error: '+str(self.global_error))
        print('Number of iterations: '+str(self.betaiter))

    def plot(self, xaxis="days"):
        if xaxis == "days":
            plt.plot(self.sim.t, self.sim.I_d, label="Simulation")
            plt.scatter(self.t_data, self.I_d_data, label="Data",color='tab:red')
            plt.xlabel("Days")
            for i in np.array(self.beta_days):
                plt.axvline(i, linestyle="dashed", color="grey")
        
        elif xaxis == "dates":
            import matplotlib.dates as mdates
            plt.plot(self.sim.dates, self.sim.I_d, label="Simulation")
            plt.scatter(self.dates_data, self.I_d_data, label="Data",color='tab:red')
            plt.xlabel("Dates")
            for i in np.array(self.beta_days):
                plt.axvline(
                    self.sim.dates[int(i)], linestyle="dashed", color="grey"
                )
            plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
            plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=10))
            plt.gcf().autofmt_xdate()
            plt.xticks(rotation=45)

        plt.ylabel("New cases")
        plt.title("Optimization results")

    def plot_step(self,step=False):
        """Plots optimization history
        """
        if not step and type(step)==bool:
            step = len(self.beta_values)
        step = min(step,len(self.beta_values))

        plt.scatter(self.t_data, self.I_d_data, label="Data",color='grey',s=20)

        if step > 0:
            beta_values = self.beta_values[:step]
            beta_days = [-np.inf] if step == 0 else self.beta_days[:step]
            beta = cv19functions.piecewise(values=beta_values,limits=beta_days) 
            sim = CV19SIM(self.cfg,I_d = self.I_d_data[0],beta=beta,mu=self.mu,inputdata=self.inputdata,**self.kwargs)
            sim.solve()
            plt.plot(sim.t,sim.I_d,label='Step: '+str(step))

            for i in beta_days:
                plt.axvline(i,linestyle='dashed',color='grey')

        plt.xlabel("Days", size=15)
        plt.ylabel("New Cases", size=15)
        plt.title('Optimization evolution')
        plt.legend(loc=0)
        plt.show()
        return

    def plot_history(self,steps=False):
        """Plots optimization history
        """
        if not steps and type(steps)==bool:
            steps = len(self.beta_values)-1
        steps = min(steps,len(self.beta_values)-1)
        
        plt.scatter(self.t_data, self.I_d_data, label="Data",color='grey',s=20)
        beta_values = [self.beta_values[0]]
        beta_days = [-np.inf]
        for i in range(steps):
            beta = cv19functions.piecewise(values=beta_values,limits=beta_days) 
            sim = CV19SIM(self.cfg,I_d = self.I_d_data[0],beta=beta,mu=self.mu,inputdata=self.inputdata,**self.kwargs)
            sim.solve()
            plt.plot(sim.t,sim.I_d,linestyle='dashed',label='Iter: '+str(i))
            beta_values.append(self.beta_values[i+1])
            if i==0:
                beta_days = [self.beta_days[i]]
                plt.axvline(self.beta_days[i],linestyle='dashed',color='grey')
            else:
                beta_days.append(self.beta_days[i])
                plt.axvline(self.beta_days[i],linestyle='dashed',color='grey')

        beta = cv19functions.piecewise(values=beta_values,limits=beta_days) 
        sim = CV19SIM(self.cfg,I_d = self.I_d_data[0],beta=beta,mu=self.mu,inputdata=self.inputdata,**self.kwargs)
        sim.solve()
        plt.plot(sim.t,sim.I_d,label='Final')

        plt.xlabel("Days", size=15)
        plt.ylabel("New Cases", size=15)
        plt.title('Optimization evolution')
        plt.legend(loc=0)
        plt.show()
        return


# Tendency Change
def tendency_change(I_d_data, t_data=None, cfg=None, sim = None, l_err_tol=100, errorfunction=LAE, startingday = 0,**kwargs):
    if not sim:
        sim = CV19SIM(config=cfg,I_d = I_d_data[0],**kwargs)
        sim.solve()
    
    # Add options for calculating local error
    l_err = errorfunction(sim.I_d,I_d_data,t_data) # Local Absolute error
    aux = np.where(np.array(l_err)>l_err_tol)[0] # first index over error tolerance
    changeindex = aux[0] if len(aux)>0 else len(I_d_data)-1

    idx = np.where(t_data>startingday)[0]
    idx = idx[0] if len(idx)>0 else len(I_d_data)-1
    return max(changeindex,idx)



# Data Management
def data_management(self, I_d_data, t_data, dates_data, tstate, initdate, cv19gmdata):
        # Data array
        if not type(I_d_data) == type(None):
            self.I_d_data = I_d_data
            self.t_data = t_data
            self.dates_data = dates_data
            self.data = None

        # State code + inital date: downloads data from remote server
        elif not type(tstate) == type(None):
            self.data = importdata(tstate=tstate, initdate=initdate)
            # Import data for seir - change this for adding more models
            print("Importing data from remote server")
            self.data.import_data_seir()
            self.I_d_data, self.t_data, self.dates_data = (
                self.data.I_d_r,
                self.data.I_d_r_tr,
                self.data.I_d_r_dates,
            )

        # cv19gmdata object
        elif not type(cv19gmdata) == type(None):
            self.data = cv19gmdata
            self.I_d_data, self.t_data, self.dates_data = (
                self.data.I_d_r,
                self.data.I_d_r_tr,
                self.data.I_d_r_dates,
            )

        # Missing data to fit
        else:
            # Raise exception
            raise Exception("Missing data to fit")

            
def beta_fit(bounds,cfg,I_d_data,t_data,beta_values=[], beta_days = [-np.inf],actualidx=-1,debug=True,global_errfunct=RMSE,**kwargs):
    starttime = time.time()
    # Find right beta for the divergence point
    bounds = [[bounds[0]],[bounds[1]]]
    opti = BETA_PIECEWISE_FIT(cfg = cfg, bounds = bounds, I_d_data = I_d_data[:actualidx],
    t_data = t_data[:actualidx],  beta_values=beta_values,beta_days=beta_days, **kwargs) # mu and inputdata are given through the kwargs

    # Choose algorithm
    algo = pg.algorithm(pg.nlopt(solver="bobyqa"))
    algo.set_verbosity(10)

    # Find parameters
    if debug:
        print("Solving for beta")
    pop = pg.population(opti, size=20)
    pop = algo.evolve(pop)

    beta_values.append(np.around(pop.champion_x[0],4))

    # Calculate Global error
    beta = cv19functions.piecewise(values=beta_values,limits=beta_days)
    
    sim = CV19SIM(cfg,I_d = I_d_data[0],beta=beta,**kwargs)
    sim.solve()

    global_errfunct = cv19errorbuild(global_errfunct)
    global_error = global_errfunct(sim.I_d, I_d_data, t_data=t_data)
    endtime = time.time()
    
    if debug:
        print('Optimization process Finished')
        print('Elapsed time:'+str(np.around(endtime - starttime,2))+' s')
        print('beta: '+str(beta_values))
        print('days: '+str(beta_days))
        print('Global Error: '+str(global_error))            
    
    return beta_values, beta_days, beta, global_error




