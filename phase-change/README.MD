# Phase change models in Basilisk solver #

The objective of this section is to develop phase-change models in Basilisk solvers in order to simulate evaporation, condensation, film boiling etc.
- The basic procedures are as follows:
    - calculating mass flux
    - solve geometric VOF equations
    - solve energy equation
    - solve N-S equation combined with mass flux in projection step
## Calculating mass flux (mass transfer models)
- Before developing mass transfer models, we use constant mass flux without heat transfer for droplet evaporation (see <a href = "https://github.com/GabrielGLK/thesis-cases/tree/master/phase-change/test_cases/mass_constant" > case-evaporation  </a> ). Because in this test case, we just focus on the numerical methods in solving VOF and N-S equation. 
    + The geometric VOF methods with mass transfer