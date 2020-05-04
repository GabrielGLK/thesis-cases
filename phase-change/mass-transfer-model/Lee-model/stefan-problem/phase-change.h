#include "tracer.h"
#include "diffusion-pc.h"
#include "curvature.h"
#include "my_function.h"

attribute {
  double D;
  double tr_eq;
  double lambda;
  double cp;
  double rho;
}

/****************** Lee-model****************************/
void lee_mass(scalar f, scalar tr, scalar m, double L_h){

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  scalar r_i[];
  foreach()
  {
    if(f[]>1e-12&&f[]<1-1e-12)
      {
        r_i[] = tr.tr_eq/(L_h*(0.5+0.5*f[])*Delta*f[]*rho1);//mass transfer intensity factor expression
        // An explicit expression of the empirical factor in a widely used phase change model
        m[] = -100*rho1*f[]*(tr[] - tr.tr_eq)/tr.tr_eq;
      }
    else
      m[] = 0;
  }
  boundary({m});
}

struct Heat_Source {
  scalar tr;
  scalar f;
  scalar m;
  double L_h;
};

mgstats heat_source (struct Heat_Source p)
{

  scalar tr = p.tr, f = p.f, m = p.m;
  double L_h = p.L_h;
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

  scalar dirichlet_source_term[];
  face vector diffusion_coef[];

  foreach()
    dirichlet_source_term[] = cm[]*m[]*L_h/(tr.rho*tr.cp);
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;

  boundary({dirichlet_source_term, diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term);
}


// diffusion the rude mass source to the neighbouring cells
struct Mass_Diffusion{
  scalar tr;
  scalar f;
  scalar m;
};

mgstats mass_diffusion (struct Mass_Diffusion p){

  scalar m = p.m, tr = p.tr, f = p.f;
  scalar mass_source[], b[], c[];
  face vector mass_diffusion_coef[];

 foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr,m});

foreach() 
  {
    mass_source[] = m[];
    b[] = -1;
    c[] = 0;
  }
boundary({mass_source, b , c});

foreach_face()
  mass_diffusion_coef.x[] = fm.x[]*sq(5*Delta);// control the mass diffusion amplitude
boundary((scalar *){mass_diffusion_coef});

return diffusion (tr, dt, D = mass_diffusion_coef, r = mass_source, beta = b, theta = c);
}