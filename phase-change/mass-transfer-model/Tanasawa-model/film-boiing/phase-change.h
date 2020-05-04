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

// delta_s function
void delta_magnini(scalar f, scalar delta ) {

  foreach()  {
    f[] = clamp(f[], 0., 1.);
    double df, df2;
    df2 = 0.;
    foreach_dimension () {
      df = (f[-1,1] + 2*f[0,1] + f[1,1] - f[-1,-1] - 2*f[0,-1] - f[1,-1])/(8*Delta);
      df2 += sq(df);
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

// Hard mass-transfer model
void mass_hard (scalar tr, scalar f, scalar m, scalar m_0, scalar m_1, double L_h)
{

  scalar delta_s[];
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

   /* step-1: initial mass flux
  */
  scalar J_evp[]; // evaporation mass flux
  /*
  the delta_s function defined here is different from div_pc[], we obtain this from 
  'CFD Modeling of Two-Phase Boiling Flows in the Slug Flow Regime with an Interface Capturing Technique' p-60
  */
  double R_int = 6.38e-8;

  // delta_s function (|\nabla c|)

  //delta_magnini(f,delta_s);
   face vector d_u[];
  // delta_s function (|\nabla c|)
  /* gradient of volume fraction*/
  foreach_face()
    //d_u.x[] = (f[-1,1] + 2*f[0,1] + f[1,1] - f[-1,-1] - 2*f[0,-1] - f[1,-1])/(2*Delta);
    d_u.x[] = (f[] - f[-1])/Delta;
  boundary((scalar *){d_u});
/* calculation of (|\nabla c|) */
  foreach()
    delta_s[] = sqrt(sq(d_u.x[]*fm.x[]) + sq(d_u.y[]*fm.x[]));
  boundary({delta_s});


  double N1 = 0, N_all = 0, N = 0;
  foreach(reduction(+:N_all) reduction(+:N1))
    {
      N1 += delta_s[]*dv();
      N_all += delta_s[]*f[]*dv();

      if (N_all > 1e-99)
        N = N1/N_all;
    }
    
  foreach()
    {
      J_evp[] = (tr[] - tr.tr_eq)/(R_int*L_h); // should be (tr[] - tr.sat)/R_int -> this is what we need to consider, not use gradient of temperature, much easier for mass flux
      m_0[] = N*J_evp[]*f[]*delta_s[];
    }
  boundary({J_evp, m_0});

/*
  step-2: smearing the mass flux 
*/

  mass_diffusion(m_1,f,m_0);

/*
  step-3: From m[] the complete source-term incorporating evap-
orative mass transfer between the phases is determined as
*/
  double N_e = 0, N_v_all = 0, N_l_all = 0, N_v = 0, N_l = 0;
  foreach(reduction(+:N_e) reduction(+:N_v_all) reduction(+:N_l_all))
  {
    N_e += m_0[]*dv();
    N_v_all += (1-f[])*m_0[]*dv();
    N_l_all += f[]*m_0[]*dv();
    if(N_v_all > 1e-99)
      N_v = N_e/N_v_all;
    if(N_l_all > 1e-99)
      N_l = N_e/N_l_all;
  }

  foreach()
  {
    if(f[] < 1e-12)
      m[] = N_v*(1-f[])*m_1[];
    else if(f[] > 1 - 1e-12)
      m[] = -N_l*f[]*m_1[];
    else
      m[] = 0;
  }
  boundary({m});

}

// Tanasawa model
//void tanasawa_model (scalar tr, scalar f, scalar m_dot, scalar m, double L_h)
void tanasawa_model (scalar tr, scalar f, scalar m_dot,  double L_h)
{
  scalar delta_s[];
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});
  delta_magnini(f,delta_s);
  double r = 0.0006;
  foreach()
  {
    if(interfacial(point,f))
      {
        m_dot[] = -(2*r/(2-r))*sqrt(20/(2*pi*8.314))*(rho2*L_h*(tr[] - tr.tr_eq))/(pow(tr.tr_eq,3/2));
        //m[] = m_dot[] *delta_s[];
      }
    else
      {
        m_dot[] = 0;
        //m[] = 0;
      }
    }
  boundary({m_dot});
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