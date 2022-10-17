#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
#import pandas as pd
import json
import logging
import numpy as np
from flask import Flask
from flask import jsonify
from flask import request
#from flask import send_file
from flask_cors import CORS

from cv19gm.cv19sim import CV19SIM
import cv19gm.utils.cv19functions as cv19functions
#import cv19gm.utils.cv19paramfit as cv19paramfit
from cv19gm.models.seir_meta import SEIRMETA
import cv19gm.utils.cv19paramfit as cv19paramfit



app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})


"""
# ------------------------- #
#     GUI communication     #
# ------------------------- #
"""
@app.route("/")
def index():
    return "<h1>Epic Suite Backend Server is up</h1>"

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
        cfg =  request.get_json(force=True)

        results = {}
        for key,value in cfg.items():
            print(key) 
            sim = CV19SIM(dict(value))
            sim.solve()
            results.update({key:sim.sims[0].results.to_json()})

        response = {'status': 'OK','results' : results}
        return jsonify(response), 200
    
    except Exception as e: 
        print(e)        
        response = {"error": "Wrong parameters"}
        return response, 400 
    
@app.route('/simulate_meta', methods=['POST'])
def simulate_meta():
    '''
    http://192.168.2.223:5003/simulate_meta?campo=1

    Estructura del código
     1.- leer parámetros
     2.- Crear objetdo de simulación
     3.- Simular. Lanzar un warning del tiempo que se puede tomar si es que es una RBM
     4.- Retornar los resultados
    '''
    try:
        cfg =  request.get_json(force=True)
        results = {}
        global_results = {}
        for key,value in cfg.items():
            print(key)
            #print(value)
            sim = SEIRMETA(dict(value))
            print('Simulating may take some time')
            sim.solve() 
            aux = {}
            for i in range(sim.nodes):
                aux[str(i)] = sim.results.loc[sim.results['node']==i].to_dict('list')               
                    
            results.update({key:json.dumps(aux)})            
            global_results.update({key:sim.global_results.to_json()})
        
        response = {'status': 'OK','results' : results,'global_results' : global_results}
        return jsonify(response), 200
    
    except Exception as e: 
        print(e)        
        response = {"error": "Wrong parameters"}
        return response, 400    

@app.route('/function', methods=['POST'])
def function():
    '''
    http://192.168.2.223:5003/function?campo=1

    Estructura del código
     1.- leer parámetros
     2.- Crear objetdo de simulación
     3.- Simular. Lanzar un warning del tiempo que se puede tomar si es que es una RBM
     4.- Retornar los resultados
    '''
    try:
        cfgfunction =  request.get_json(force=True)
        function = cv19functions.build(cfgfunction['function'])
        t = np.linspace(cfgfunction['t_init'],cfgfunction['t_end'],(cfgfunction['t_end']-cfgfunction['t_init'])*10+1)
        functionarray = function(t)
        response = {'status': 'OK','results' : {'t':list(t),'function':list(functionarray)}}
        return jsonify(response), 200
    
    except Exception as e: 
        print(e)        
        response = {"error": "Wrong parameters"}
        return response, 400    


@app.route('/datafit', methods=['POST'])
def datafit():
    """Optimal parameters for fitting data

    Returns:
        _type_: _description_
    """
    #if "cv19gm.utils.cv19paramfit" not in sys.modules:
        
        #print('You have not imported the {} module'.format(modulename))
    
    input =  request.get_json(force=True)

    results = {}
    
    fit = cv19paramfit.SEQUENTIAL_FIT(cfg = 'SEIR.toml', I_d_data=np.array(list(json.loads(input['I_d_data']).values())),t_data = np.array(list(json.loads(input['t_data']).values())),global_errortol=200, local_errortol=250, 
        intervalsize=10, maxintervals=5, bounds_beta=[0,1], bounds_mu=[0,4],tE_I = input['tE_I'],tI_R = input['tI_R'])
    
    fit.optimize()
        
    
    results['beta_values'] = str(fit.beta_values)
    results['beta_days'] = str(fit.beta_days)
    results['mu'] = fit.mu
    results['simulation'] =  fit.sim.sims[0].results.to_json()
    
    response = {'status': 'OK','results' : results}
    try:  
        return jsonify(response), 200
    
    except Exception as e: 
        print(e)        
        response = {"error": "Wrong parameters"}
        return response, 400    


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--debug", action='store_true', help = "debug mode")
    args = parser.parse_args()
    if args.debug:
        debug = True
    else:
        debug = False

    if debug:
        print("Initializing server in debugging mode")
        app.run(host="0.0.0.0", port=5003, debug=True)
    else: 
        print("Initializing server in production mode")
        from waitress import serve
        serve(app, host="0.0.0.0", port=5003)
    
    
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



