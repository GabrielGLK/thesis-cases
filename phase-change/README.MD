# Phase change models in Basilisk solver #

The objective of this section is to develop phase-change models in Basilisk solvers in order to simulate evaporation, condensation, film boiling etc.
- The basic procedures are as follows:
    - calculating mass flux
    - solve geometric VOF equations
    - solve energy equation
    - solve N-S equation combined with mass flux in projection step
## Calculating mass flux (mass transfer models)
