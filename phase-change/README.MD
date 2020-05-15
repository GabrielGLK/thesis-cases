# Phase change models in Basilisk solver #

The objective of this section is to develop phase-change models in Basilisk solvers in order to simulate evaporation, condensation, film boiling etc.
- The basic procedures are as follows:
    - calculating mass flux
    - solve geometric VOF equations
    - solve energy equation
    - solve N-S equation combined with mass flux in projection step
## Calculating mass flux (mass transfer models)
- Before developing mass transfer models, we use constant mass flux without heat transfer for droplet evaporation (see <a href = "https://github.com/GabrielGLK/thesis-cases/tree/master/phase-change/test_cases/mass_constant" > case-evaporation  </a> ). Because in this test case, we just focus on the numerical methods in solving VOF and N-S equation. 
    - The geometric VOF methods with mass transfer
        In Basilisk, the geometric VOF method is adopt to track two immicible phase interface, which including interface reconstruction (at beginning of each time step) and interface advection (consistent VOF flux and WY scheme is used). One thing needs to be noted that in interface advection step, the prediction of fluid is divergence-free. Therefore, the VOF methods maybe not available directly to the existed VOF method due to the phase change source term. Actually, people more like to use level-set method to simulate phase-cahnge problems, and the divergence-free problem in VOF equation sometime not considered. Here, we use the idea of Leon Malan to extend the liquid domain to vapor domain and the difference between origin vapor velocity and ghost domain velocity would produce a divergene-free velocity for the entire computational domain.
        ![equation](https://latex.codecogs.com/gif.latex?\dpi{400}\frac{\partial c}{\partial t}   \approx \frac{c^{n+1}-c^n}{\Delta t} = \frac{1}{V}\int_V \frac{\partial H_1}{\partial t} + \frac{\partial H_2}{\partial t} dV)