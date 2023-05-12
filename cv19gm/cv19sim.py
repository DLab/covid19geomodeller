import toml

"""
CV19SIM interphaces the user with all the cv19gm library tools.
"""

def CV19SIM(config = None, compartmentalmodel=None, verbose=False,**kwargs):
        
        # Read the configuration file
        if not config:
            if not compartmentalmodel:
                raise Exception("Missing compartmental model definition")
            if verbose:
                print("Using default configuration file for "+compartmentalmodel+" model")            
        else:
            # TOML file
            if not type(config) == dict:
                compartmentalmodel = toml.load(config)['model']['name']
            # Dictionary
            else:
                compartmentalmodel = config['model']['name']
           
        
        if compartmentalmodel == 'SEIR':            
            from cv19gm.models.seir import SEIR
            return SEIR(config,verbose=verbose,**kwargs)
             
        elif compartmentalmodel == 'SEIRHVD':            
            from cv19gm.models.seirhvd import SEIRHVD
            return SEIRHVD(config,verbose=verbose,**kwargs)

        elif compartmentalmodel == 'SIR':            
            from cv19gm.models.sir import SIR
            return SIR(config,verbose=verbose,**kwargs)
            
        elif compartmentalmodel == 'SEIRTQ':
            from cv19gm.models.seirtq import SEIRTQ
            return SEIRTQ(config,verbose=verbose,**kwargs)
                        
        elif compartmentalmodel == 'SEIR_META':                
            from cv19gm.models.seir_meta import SEIRMETA
            return SEIRMETA(config,verbose=verbose,**kwargs)
            
        else:
            raise('Incorrect model: '+str(compartmentalmodel))
