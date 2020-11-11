# Examination Dynamics Study

This work studies how different examinations campings can affect a pandemic evolution. In order to facilitate its lecture, the work is split in different jupyter notebooks

## Model Definition

\begin{equation} 
\frac{dS}{dt} = -\beta \alpha \frac{S I}{N+k_{I}I+k_{R}R}
\end{equation}

\begin{equation} 
\frac{dE}{dt} = \beta \alpha \frac{S I}{N+k_{I}I+k_{R}R} - \sigma E  
\end{equation}
\begin{equation} 
\frac{dI}{dt} = \sigma E - \gamma I - \lambda_{e} \sqcap(\omega_{e} t) P_{I}
\end{equation}
\begin{equation} 
\frac{dR}{dt} = \gamma I + \lambda_{e} \sqcap(\omega_{e} t) P_{I}
\end{equation}
\begin{equation} 
S_i+E_i+I_i+R_i = N_i
\end{equation}
\begin{equation} 
P_{I} = \frac{I}{N}
\end{equation}
\begin{equation} 
e(t) = \lambda_{e} \sqcap(\omega_{e} t)
\end{equation}_
\begin{equation} 
e_{I}(t) = \lambda_{e} \sqcap(\omega_{e} t)
\end{equation}

## Epidemiological Parameters
* **beta:** Infection rate

* **mu:** Initial exposed obtained from the initial infected mu=E0/I0

* **Scale Factor:** Proportion of real infected compared to reported ones (1: all the infecteds are reported)

* **Sero Prevalence Factor:** Adjust the proportion of the population that enters the virus dynamics

* **Exposed Infection:** rate compared to the infected (0 the don't infect, 1 the infect in the same rate as the infected )

  

## Table of Contents:
* Single Simulation Study
    * Parameter Settings
    * SEIR plot
    * Accumulated Infected
    * New Daily Infected
    * Ammount of exams
* Sensitivity Analysis: 
    * Examination Rate
        * Peak Size
        * Peak time
        * Prevalence
    * Examination Period
        * Peak Size
        * Peak time
        * Prevalence
    * Examination Duty
        * Peak Size
        * Peak time
        * Prevalence  
    * Examination Accuracy
        * Peak Size
        * Peak time
        * Prevalence        
* Multidimensional Analysis:
    * Examination rate vs Examination periods
    * Examination Rate vs Duty