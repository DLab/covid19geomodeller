#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

from timeit import default_timer as timer

from datetime import datetime

import logging
import json

from cv19gm.cv19sim import CV19SIM

from flask import Flask
from flask import jsonify
from flask import request
from flask import send_file
from flask_cors import CORS

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})


"""
# ------------------------- #
#     GUI communication     #
# ------------------------- #
"""

@app.route('/health_check', methods=['GET'])
def health_check():
    '''
    health_check
    '''
    app.logger.info("health_check")
    response = {'status': 'OK'}
    return jsonify(response), 200

@app.route('/test', methods=['POST'])
def test():
    try: 

        # 1. Get params
        form =  request.form
        data = form.to_dict()
        print(data)
        print(type(data))

        # 2. Build cv19sim object
        sim = CV19SIM(data)

        # 3. Simulate
        sim.solve()

        response = {'status': 'OK','result' : data}

        return jsonify(response), 200
    except:
        response = {"error": "Unkown error"}
        return response, 200

@app.route('/simulate', methods=['POST'])
def simulate():
    '''
    http://192.168.2.223:5003/simulate?campo=1

    Estructura del código
     1.- leer parámetros
     2.- Crear objetdo de simulación
     3.- Simular. Lanzar un warning del tiempo que se puede tomar si es que es una RBM
     4.- Retornar los resultados
    '''
    try: 

        # 1. Get params
        form =  request.form
        cfg = form.to_dict()
        #print(data)
        #print(type(data))


        # 2. Build cv19sim object
        sim = CV19SIM(cfg)

        # 3. Simulate
        sim.solve()

        response = {'status': 'OK','result' : sim.results.to_json()}

        return jsonify(response), 200
    except:
        response = {"error": "Unkown error"}
        return response, 200    

@app.route('/refineuni', methods=['POST'])
def datafit():
    '''
    # ----------------------------------------------------------- #
    #         Find optimal parameters for fitting the data        #
    # ----------------------------------------------------------- #
    '''
    try: 
        print("Optimizing parameteres")


        response = "test"
        #print(response)
        return jsonify(response), 200

    except Exception as e: 
        print(e)
        response = {"error": str(e)}
        return response, 200


@app.route('/SEIR2refine', methods=['GET'])
def SEIR2refine():
    '''
    # ----------------------------------------------------------- #
    #           Parameters Refine for intercomunal model          #
    # ----------------------------------------------------------- #
    '''
    print("Refine")
    try: 
        # Manage input parameters
        if request.form:
            print("I have form data")
            state = request.form.get('state') #State (Region)
            comuna = request.form.get('comuna') # District 
            qp = int(request.form.get('qp')) # Quarantine period
            mov = float(request.form.get('mov')) # Quarantine movilty
            tsim = int(request.form.get('tSim')) # Simulation time
            movfunct = str(request.form.get('movfunct'))    

        if request.json:
            print("I have json")        
            state = request.json['state']
            comuna = request.json['comuna']
            qp = int(request.json['qp'])
            mov = float(request.json['mov']) # Quarantine movilty
            tsim = int(request.json['tSim'])
            movfunct = str(request.json['movfunct'])
        
        #tsim = 100
        #mov = 0.2          
        if qp == -1:
            qp = tsim
             
        # ---------------------------- #
        #        Get Refine Data       #
        # ---------------------------- #
        path = '../Data/unirefine/'
        
        # Get data for different quarantine periods
        # Import data, parameters and initdate

        #parameters = pd.read_csv(path+'parameters_qp0.csv',index_col=0)
        #init_date = pd.read_csv(path+'initdate_qp0.csv',index_col=0)
        parameters = pd.read_csv(path+'parameters.csv',index_col=0)
        init_date = pd.read_csv(path+'initdate.csv',index_col=0)        

        ## Find data for paramters given
        parameters = parameters[comuna]
        init_date = init_date[comuna][0]
        results = SDSEIR.simulate(state,comuna,parameters[0],parameters[1],parameters[2],parameters[3],qp = qp,mov = mov,tsim = tsim, movfunct=movfunct)
        S = results['S'].tolist()
        E = results['E'].tolist()
        I = results['I'].tolist()
        R = results['R'].tolist()
        t = list(results['t'])
        
        Ti = 1/parameters[1]
        Tr = 1/parameters[2]

        # Round results:
        S = [round(i) for i in S]
        E = [round(i) for i in E]
        I = [round(i) for i in I]
        R = [round(i) for i in R]

        #print(results)
        print("done simulating")    

        response = {'status': 'OK','S':S,'E':E,'I':I,'R':R,'t':t, 'Ti':Ti,'Tr':Tr,'beta':parameters[0],'r0':parameters[3],'init_date':init_date,'I_peak':max(I),'R_total':(max(R))}
        #print(response)
        return jsonify(response), 200

    except Exception as e: 
        print(e)
        response = {"error": str(e)}
        return response, 200


@app.route('/SEIR2simulate', methods=['GET'])
def SEIR2simulate():
    '''
    # ----------------------------------------------------------- #
    #             Multidistricts SEIR Model simulation            #
    # ----------------------------------------------------------- #
    '''
    
    # ------------------------ #
    #      Input parameters    #
    # ------------------------ #
    try: 
        if request.form:
            print("I have form data")
            state = str(request.form.get('state'))
            urbancenter = str(request.form.get('urbancenter'))
            qp = int(request.form.get('qp'))
            beta = float(request.form.get('beta'))
            Ti = int(request.form.get('Ti')) #sigma-1
            Tr = int(request.form.get('Tr')) #gamma-1
            r0 = float(request.form.get('r0'))    
            mov = float(request.form.get('mov')) #Quarantine movilty
            tsim = int(request.form.get('tSim')) 
            movfunct = str(request.form.get('movfunct'))                   

        if request.json:
            print("I have json")        
            state = str(request.json['state'])
            urbancenter = str(request.json['urbancenter'])
            qp = int(request.json['qp'])
            beta = float(request.json['beta'])
            Ti = int(request.json['Ti']) #sigma-1
            Tr = int(request.json['Tr']) #gamma-1
            r0 = float(request.json['r0'])
            mov = float(request.json['mov'])
            tsim = int(request.json['tSim'])
            movfunct = str(request.json['movfunct'])

        # Get Comunas from OD 
        # falta identificar centros urbanos     
        comunas = []
        if qp ==-1:
            qp = tsim

        if Ti==0:
            sigma = 1
        else:
            sigma = 1/Ti

        if Tr==0:
            gamma = 1
        else:
            gamma = 1/Tr    
        tci = None
        #params = beta sigma gamma mu
        params = [beta,sigma,gamma,r0] 
        
        # inputs: state, centrourbano, beta, sigma, gamma, mu,qp, mov,tsim, tci, movfunct 
        # outputs: {'S':sim.S,'E':sim.E,'I':sim.I,'R':sim.R,'tout':sim.t}
        results = MDSEIR.simulate_multi(state,comunas,params,qp,mov,tsim,tci,movfunct)
        S = results['S'].tolist()
        E = results['E'].tolist()
        I = results['I'].tolist()
        R = results['R'].tolist()
        t = list(results['tout'])

        # Round results:
        S = [round(i) for i in S]
        E = [round(i) for i in E]
        I = [round(i) for i in I]
        R = [round(i) for i in R]
                
        init_date = results['init_date']
        response = {'status': 'OK','S':S,'E':E,'I':I,'R':R,'t':t,'init_date':init_date,'I_peak':max(I),'R_total':(max(R))}
        return jsonify(response), 200

    # Solve
    #integr(self,t0,T,h,E0init=False)
    except Exception as e:
        print(e)
        response = {'status': 'OK',"error": str(e)}
        return response, 200



if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5003, debug=True)
else:
    # setup logging using gunicorn logger
    formatter = logging.Formatter(
        '[%(asctime)s.%(msecs)03d] [%(name)s] [%(levelname)s] - %(message)s',
        '%d-%m-%Y %H:%M:%S'
    )
    gunicorn_logger = logging.getLogger('gunicorn.error')
    app.logger.handlers = gunicorn_logger.handlers
    app.logger.handlers[0].setFormatter(formatter)
    app.logger.setLevel(logging.DEBUG)



