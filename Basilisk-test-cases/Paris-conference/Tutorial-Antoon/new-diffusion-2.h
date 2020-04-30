/*
make a function to calculate the fluxes(F) of c in a face vector field.
We also need to provide the difusivity on faces, which could be constant
*/
#include "poisson.h"
/*
struct Diffusion {
  // mandatory
  scalar f;
  double dt;
  // optional
  face vector D;  // default 1
  scalar r, beta; // default 0
  scalar theta;   // default 1
};

trace
mgstats diffusion (struct Diffusion p)
{
    if (p.dt == 0.) {
    mgstats s = {0};
    return s;
  }

   scalar f = p.f, r = automatic (p.r);

   const scalar idt[] = - 1./p.dt;
  (const) scalar theta_idt = p.theta.i ? p.theta : idt;
  
  if (p.theta.i) {
    scalar theta_idt = p.theta;
    foreach()
      theta_idt[] *= idt[];
  }

  if (p.r.i)
    foreach()
      r[] = theta_idt[]*f[] - r[];
  else // r was not passed by the user
    foreach()
      r[] = theta_idt[]*f[];

 
  scalar lambda = theta_idt;
  if (p.beta.i) {
    scalar beta = p.beta;
    foreach()
      beta[] += theta_idt[];
    lambda = beta;
  }
  boundary ({lambda});

 return poisson (f, r, p.D, lambda);
}      
*/


void flux_diffusion(scalar c, (const) face vector D, face vector F){
    boundary({c});

    foreach_face()
        F.x[] = -D.x[]*(c[] - c[-1])/Delta;
}

/*
Using this function, we can easily compute the tendency for c in a cell(j) with N faces
*/

void tendency_from_flux(face vector F, scalar dc){
    boundary_flux({F}); // flux on levels
    foreach(){
        dc[] = 0;
        foreach_dimension()
            dc[] += (F.x[] - F.x[1])/Delta; // similar to VOF method
    }
}
// advance in time step
void advance(scalar c, scalar dc, double dt){
    foreach()
        c[] += dt*dc[];
}



struct DIFFUSION{
    scalar c;
    double dt;
    face vector D;
};

trace
mgstats diffusion(struct DIFFUSION p){
    scalar c = p.c;
    face vector D = p.D;
    double dt = p.dt;

    face vector F[];
    scalar dc[];
    flux_diffusion(c,D,F);
    tendency_from_flux(F,dc);
    advance(c,dc,dt);
}



/*
// now it's time to define our own diffusion solver
void diffusion(scalar c, double dt, face vector D){
    face vector F[];
    scalar dc[];
    flux_diffusion(c,D,F);
    tendency_from_flux(F,dc);
    advance(c,dc,dt);
}
*/

// second-order accurate in time
void diffusion_second (scalar c, double dt, face vector D){
    face vector F[];
    scalar dc[], c_temp[];
    foreach()
        c_temp[] = c[];
    flux_diffusion(c_temp,D,F);
    tendency_from_flux(F,dc);
    advance(c_temp,dc,dt/2);

    flux_diffusion(c_temp,D,F);
    tendency_from_flux(F,dc);
    advance(c,dc,dt);

}
