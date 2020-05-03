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

// constant mass source scalar function
void mass_source(scalar f, scalar m)
{
  foreach()
  {
    if(f[]>1e-12&&f[]<1 - 1e-12)
      m[] = 0.05;
    else
      m[] = 0;
  }
  boundary({f,m});
}

// delta_s function
void delta_magnini(scalar f, scalar delta ) {

  foreach()  {
    f[] = clamp(f[], 0., 1.);
    double df, df2;
    df2 = 0.;
    df = 0;
    foreach_dimension () {
      df += (f[-1,1] + 2*f[0,1] + f[1,1] - f[-1,-1] - 2*f[0,-1] - f[1,-1])/(8*Delta);
      df2 = sq(df);
    }
    delta[] = sqrt(df2);
  }
  boundary({delta});
}

void smoother(scalar fr, scalar ff) {
  foreach()  
    ff[] = (fr[-1,1] + 2*fr[-1,0] + fr[-1,-1] + fr[1,1] + 2*fr[1,0] + fr[1,-1] + 2*fr[0,1] + 4*fr[0,0] + 2*fr[0,-1])/16;
  boundary({fr,ff});
}



/****************** Lee-model****************************/
void lee_mass(scalar f, scalar m, scalar tr){

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

  foreach(){
    if(f[]>1e-12&&f[]<1-1e-12)
      m[] = 5000*f[]*(tr[] - tr.tr_eq)/(tr.tr_eq+373);
    else
      m[] = 0;
  }
  boundary({m});
}

/****************** mass transfer model-4 (calculate phase-change velocity -> mass transfer rate)****************************/
void mass_dot_4 (scalar f, scalar tr,  scalar m, scalar mm, double L_h) {
  face vector v_pc[];

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  
  face vector gtr[];
  foreach_face()
    gtr.x[] = (tr[] - tr[-1])/Delta;
  boundary((scalar*){gtr});

  vector n[];
  compute_normal (f, n);
  
  foreach_face() {
    v_pc.x[] = 0.;

    if (interfacial(point, f) || interfacial(neighborp(-1), f) ) {
      coord nf;
      foreach_dimension()
        nf.x = 0.;
      if (interfacial(point, f)) {
        foreach_dimension()
          nf.x += n.x[];
      }
      if (interfacial(neighborp(-1), f)) {
        nf.x += n.x[-1];
        nf.y += n.y[-1];
      }
   

      double norm = 0.;
      foreach_dimension()
        norm += fabs(nf.x);
      foreach_dimension()
        nf.x /= (tr.inverse ? -norm : norm);
      
      
      if (nf.x > 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1, 0] + fabs(nf.y)*(nf.y > 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[-1, 0]+ fabs(nf.y)*(nf.y > 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  }
  boundary((scalar *){v_pc});
  foreach()
    {
      m[] = 0;
      mm[] = 0;
    foreach_dimension()
    {if(interfacial(point,f))
      {
        m[] += tr.lambda*v_pc.x[]/L_h*fabs((f[]-f[-1])/Delta);
        mm[] += tr.lambda*v_pc.x[]/L_h*n.x[];
      }
     else
      {
        m[] += 0;
        mm[] += 0;
      }
    }
    }
  boundary({m});
}

/****************** equilibrium saturated temperature on interface *************************/
struct Dirichlet_Diffusion {
  // mandatory
  scalar tr;
  scalar f;
  int max_level;
  double dt;
  double time_factor;
  // optional
  scalar tr_op; // default uniform 1.
};

mgstats dirichlet_diffusion (struct Dirichlet_Diffusion p) {
  
  scalar tr = p.tr, f = p.f, tr_op = automatic (p.tr_op);
  int max_level = p.max_level;
  double time_factor = p.time_factor;

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  if (p.tr_op.i)
    boundary ({tr_op});

  scalar volumic_metric[], dirichlet_source_term[], dirichlet_feedback_term[];
  face vector diffusion_coef[];

  foreach() {
    volumic_metric[] = cm[];
    if (p.tr_op.i)  
      dirichlet_source_term[] = cm[]*tr_op[]*tr.tr_eq*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
    else
      dirichlet_source_term[] = cm[]*tr.tr_eq*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
    dirichlet_feedback_term[] = - cm[]*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
  }
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;

  boundary({volumic_metric, dirichlet_source_term, dirichlet_feedback_term, diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term,
			              beta = dirichlet_feedback_term, theta = volumic_metric);
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
  mass_diffusion_coef.x[] = fm.x[]*sq(1.5*Delta);// control the mass diffusion amplitude
boundary((scalar *){mass_diffusion_coef});

return diffusion (tr, dt, D = mass_diffusion_coef, r = mass_source, beta = b, theta = c);
}