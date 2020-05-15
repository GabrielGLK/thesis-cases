#include "tracer-pc.h"
#include "diffusion-pc.h"
#include "curvature.h"
#include "my_function.h"
#include "embed-pc.h"

attribute {
  double D;
  double tr_eq;
  double lambda;
  double cp;
  double rho;
  scalar m;
}

// ghost cell method
double interpolate_1 (Point point, scalar s, coord p)
{
  int i = sign(p.x), j = sign(p.y);
  if (f[i] && f[0,j] && f[i,j])
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (f[i])
	val += fabs(p.x)*(s[i] - s[]);
      else if (f[-i])
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
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


/***********************************************************************/
/***********************************************************************/
/************************mass transfer models***************************/
/****************** Lee-model****************************/
void lee_mass(scalar tr, scalar f, scalar m, double L_h){

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  scalar r_i[];
  foreach()
  {
    coord n = normal(point,f);
    if(f[]>1e-12&&f[]<1-1e-12)
      {
      r_i[] = cm[]*tr.lambda*tr.tr_eq/(L_h*(0.5+0.5*f[])*Delta*(n.x + n.y)*f[]*rho1);//mass transfer intensity factor expression
      // An explicit expression of the empirical factor in a widely used phase change model
      m[] = -r_i[]*rho1*f[]*(tr[] - tr.tr_eq)/tr.tr_eq;
      }
    else
      m[] = 0;
  }
  boundary({m});
}


/****************** Tanasawa-model****************************/
//void tanasawa_model (scalar tr, scalar f, scalar m_dot, scalar m, double L_h)
void tanasawa_model (scalar tr, scalar f, scalar m_dot,  double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});
  double r = 0.5;
  foreach()
  {
    if(f[]>1e-12&&f[]<1-1e-12)
      {
        m_dot[] = cm[]*(2*r/(2-r))*sqrt(18/(2*pi*8.314))*(rho2*L_h*(tr[] - tr.tr_eq))/(pow(tr.tr_eq,3/2));
      }
    else
      {
        m_dot[] = 0;
      }
    }
  boundary({m_dot});
}

/*************************  Hard mass-transfer model ***********************/
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

/****************** sharp interface-model****************************/
void sharp_simple_model_2 (scalar tr, scalar f, scalar temp,  face vector v_pc, double L_h) {

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

    if ((interfacial(point, f) || interfacial(neighborp(-1), f))) {
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
        nf.x /= norm;
      
      
      if (nf.x > 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1, 0] + fabs(nf.y)*(nf.y > 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[-1, 0]+ fabs(nf.y)*(nf.y > 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  } 
  boundary((scalar *){v_pc});
  double h = 0;
  face vector fs[];
  face_fraction(f,fs);
  scalar t_sat[];
  scalar d[], tep[];
  foreach()
    { temp[] = 0;
      d[] = 1-f[];
      if(interfacial(point,f)){
        t_sat[] = tr.tr_eq;
        double df = 0, tep = t_sat[];
        coord n  = facet_normal( point, f ,fs) , p;
        normalize(&n);
        double alpha  = plane_alpha (f[], n);
        line_length_center (n, alpha, &p);
        h = dirichlet_gradient (point, tr, f, n, p, tep, &df) ;
        temp[] = tr.lambda*h/L_h;
    }
    }
  boundary({temp});
}

void sharp_simple_model (scalar tr, scalar f, scalar temp,  face vector v_pc, double L_h) {

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

    if ((interfacial(point, f) || interfacial(neighborp(-1), f))) {
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
        nf.x /= norm;
      
      
      if (nf.x > 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1, 0] + fabs(nf.y)*(nf.y > 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[-1, 0]+ fabs(nf.y)*(nf.y > 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
    }
  } 
  boundary((scalar *){v_pc});
  scalar dd[];
  foreach()
    { 
      if(interfacial(point,f))
        {
          foreach_dimension()
            temp[] += tr.lambda*v_pc.x[]*n.x[]/L_h;
            }
      else
          temp[] = 0;
    }
  boundary({temp});
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
    dirichlet_source_term[] = cm[]*m[]*L_h/(tr.cp*tr.rho);
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;

  boundary({dirichlet_source_term, diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term);
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