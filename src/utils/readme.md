# Data Libraries (Work in progress)
## Variable Names construction:
Variable names are constructed as follows:
**Var_state_real_time/add** 
where:
* Var: Variable name
* state: if it's active **(act)**, accumulated **(acc)** or daily **(d)**
* real: if it's real data it will be explicited with an **r**. Simulated data doesn't show anything
* time: if the variable is real data it has 2 time vectors:
  * **dates**: Real time data with a datetime.datetime objects
  * **tr**: Relative day since the beginning of the simulation
 * add: Additional information for some variables. For example here we can differentiate between suspected and confirmed deaths as follow. D_acc_r_suspected  if the variable is real data it has 2 time vectors:
  * **dates**: Real time data with a datetime.datetime objects
  * **tr**: Relative day since the beginning of the simulation
  
### Examples:
Infected accumulated data:
* I_acc: Simulation data
* I_acc_r: Real data
* I_acc_r_dates: datetime data
* I_acc_r_tr: relative days

* **Ir**
* **I_d**
* **I_ac**
* **Br**
## Data Sources:
### Infected

### Deaths
#### Published by Minsal

#### Published by DEIS

### Hospitalized

### Lockdowns and efemerides
#### Lockdowns

#### Efemerides: 
Important dates that might affect pandemic evolution 

### Exams



### Geography and Demography
* Population
* Adjacency 
