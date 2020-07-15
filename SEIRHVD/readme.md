# Plataforma de simulación de modelos compartimentales

asdfasd

# Instalación

Github

Instalación de librerías

python

vpn



# Modelo SEIRHVD



## Parámetros Generales


## Parámetros del Modelo
### Región Por CUT 
tstate = '13'
#### Fecha Inicial
initdate = datetime(2020,5,15)

### Parametros EPI
beta: tasa de contagio = 0.117 
mu: Razon E0/I0
ScaleFactor: Factor de Escala, Numero de infectados por sobre los reportados
SeroPrevFactor: Factor de Seroprevalencia. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection: Proporcion en la que contagian los expuestos
tsim: Tiempo de simulacion

## Creación de Escenarios

inputarray = np.array([tsim,max_mov,rem_mov,qp,iqt,fqt,movfunct])

* tsim: Tiempo de simulación
* max_mov: Movilidad máxima durante tiempo sin cuarentena
* rem_mov: Movilidad remanente durante tiempo de cuarentena
* qp: Periodo de Cuarentena para cuarentenas alternantes - qp días con cuarentena, luego qp días sin cuarentena
* iqt: Día de inicio de cuarentena (desde el inicio de la simulación)
* fqt: Día de fin de cuarentena (desde el inicio de la simulación)
* movfunct: Función de movilidad
    * 0: Cuarentena total durante el período comprendido entre iqt y fqt
    * 1: Cuarentena alternante tipo onda Cuadrada con período qp
    * 2: Cuarnetena tipo diente de cierra con período qp






# Glosario de Variables:
    # initial params
    initdate = None
    tsim = None
    May15 = None
    tstate = None
    
    # Parameters
    beta = None
    mu = None
    ScaleFactor = None
    SeroPrevFactor = None
    expinfection = None
    
    # Quarantine Scenarios
    inputarray  = None 
    numescenarios  = None 
    
    # Data
    # Sochimi
    sochimi  = None 
    Hr  = None 
    Vr  = None  
    Vr_tot  = None 
    Hr_tot  = None  
    sochimi_dates  = None 
    sochimi_tr  = None 
    # Infected minsal
    Ir  = None 
    Ir_dates  = None 
    tr  = None 
    # Death minsal
    Br  = None 
    Br_dates  = None 
    Br_tr  = None 
    # Death excess
    ED_RM  = None 
    ED_RM_dates  = None 
    ED_tr  = None 
    ED_RM_ac  = None 
    # Death deis 


    # Sim and Auxvar:
     sims = None  
    Vcmodel = None  
    Hcmodel = None  
    sims = None  
    Vmax = None  
    Hmax = None  
    T  = None 
    S  = None
    H_sum  = None
    H_bed  = None
    H_vent  = None
    Iac  = None
    I  = None
    I_act  = None
    I_as  = None
    I_mi  = None
    I_se  = None
    I_cr  = None
    I_sum  = None
    E  = None
    E_as  = None
    E_sy  = None
    B  = None
    D  = None
    R  = None
    V  = None
    t  = None
    dt  = None
    idx  = None
    H_crin  = None
    H_in  = None
    H_out  = None
    H_sum  = None
    H_tot  = None
    CH  = None
    CV  = None
    ACH  = None
    ACV  = None
    peakindex  = None
    peak  = None
    peak_t  = None
    peak_date  = None
    population  = None
    infectedsusc  = None
    infectedpop  = None
    err_bed  = None
    err_vent  = None
    err_Iactives  = None
    H_colapsedate  = None
    V_colapsedate  = None


# Funciones de ploteo

## Función Standard

simulation.pĺotvariable(enddate =  datetime(2020,7,30),days=0, reales= True,ylim = 0,norm=1,scalefactor = False)



## Lista de funciones: 

### plotactivos()



# Creación de Tablas 





## Parámetros de ajuste

## Santiago

### Ajuste por ventiladores

beta = 0.117 #0.25#0.19 0.135
mu = 0.6 #2.6 0.6
ScaleFactor = 1.9 #4.8
SeroPrevFactor = 0.5 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos



### Ajuste por Muertos minciencia (SEIRHVD Class 2)

beta = 0.117 #0.25#0.19 0.135
mu = 0.15 #2.6 0.6
ScaleFactor = 2.25 #4.8
SeroPrevFactor = 0.5 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos



### Ajuste por Muertos minciencia para adelantar peak(SEIRHVD Class 2)

beta = 0.115 #0.25#0.19 0.135
mu = 0.15 #2.6 0.6
ScaleFactor = 4 #4.8
SeroPrevFactor = 0.3 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos



### Ajuste por exceso de Muertos (SEIRHVD Class 2)

beta = 0.14 #0.25#0.19 0.135
mu = 0.01 #2.6 0.6
ScaleFactor = 2.0 #4.8
SeroPrevFactor = 0.5 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos



beta = 0.119 #0.25#0.19 0.135
mu = 0.7#2.6 0.6
ScaleFactor = 2.5 #4.8
SeroPrevFactor = 0.5 # Sero Prevalence Factor. Permite ajustar la cantidad de gente que entra en la dinamica
expinfection = 1 # Proporcion en la que contagian los expuestos


### Ajuste por datos de Camas

#### Optimos fit de camas 13-04

beta = 0.117
mu = 0.6
fI = 1.9
SeroPrevFactor = 0.5

## Atacama




# 