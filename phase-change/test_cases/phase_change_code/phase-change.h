#include "tracer-pc.h"
#include "diffusion-pc.h"
#include "curvature.h"
#include "my_function.h"

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
      m[] = 0.5;
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
    mass_source[] = cm[]*m[];
    b[] = -1;
    c[] = 0;
  }
boundary({mass_source, b , c});

foreach_face()
  mass_diffusion_coef.x[] = fm.x[]*sq(2*Delta);// control the mass diffusion amplitude
boundary((scalar *){mass_diffusion_coef});

return diffusion (tr, dt, D = mass_diffusion_coef, r = mass_source, beta = b, theta = c);
}

struct Volume_Diffusion{
  scalar tr;
  scalar f;
};

mgstats volume_diffusion (struct Volume_Diffusion p){

  scalar tr = p.tr, f = p.f;
  scalar mass_source[], b[], c[];
  face vector mass_diffusion_coef[];

 foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});

foreach() 
  {
    mass_source[] = f[];
    b[] = -1;
    c[] = 0;
  }
boundary({mass_source, b , c});

foreach_face()
  mass_diffusion_coef.x[] = fm.x[]*sq(1*Delta);// control the mass diffusion amplitude
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
    r_i[] = 27.2*tr.tr_eq/(L_h*(0.5+0.5*f[])*Delta*f[]*rho1);//mass transfer intensity factor expression

    if(interfacial(point,f))
      m[] = 1*rho1*(f[] + f[1])/2*(tr[] - tr.tr_eq)/tr.tr_eq;
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
  //double a = 0.2;// because we don't know the liquid molecular weight, so we give the guess of mass transfer intensity multipy liquid molecular weight
  double b = 0.01;
  foreach()
  {
    if(interfacial(point,f))
      {
        if(f[] != f[1])
          m_dot[] = cm[]*(2*b/(2-b))*sqrt(18/(2*pi*8.314))*(rho2*L_h*(tr[] - tr.tr_eq))/(pow(tr.tr_eq,3/2));
        else
          m_dot[] = 0;
      }
    }
  boundary({m_dot});
}

/****************** zhang-model****************************/
// This mass transfer model from Jie zhang's paper, calculating mass flux
void zhang_model (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  scalar f_v[];
  foreach()
    {
      m[] = 0;
      f_v[] = 1 - f[]; // change origin volume fraction in order to calculate temperature variation in vapor cell
      coord n  = interface_normal( point, f) , q, s; // interface normal, vapor interface normal, interface center coordinate
      q.x = -n.x, q.y = -n.y; // vapor interface normal direction
      double alpha_new  = plane_alpha (f_v[], q); // distance in vapor portion in interfacial cell
      double alpha  = plane_alpha (f[], n); // distance in liquid portion in interfacial cell
      plane_area_center(q, alpha_new, &s);// interface center
      double xc = 0, yc = 0;
      double delta_d = 1;// according to the paper, this value is between 1 and 2
      if(interfacial(point,f))
        {
          xc += s.x*Delta-n.x*delta_d*Delta;
          yc += s.y*Delta-n.y*delta_d*Delta;
          coord pp;
          pp.x = xc, pp.y = yc;// mirror point in local corner neighboring cell
          double d1 = fabs(n.x*xc + n.y*yc - alpha*Delta - Delta*0.5*(n.x+n.y))/sqrt(sq(n.x) + sq(n.y));// distance for temperature gradient
          double tt = interpolate_1 (point, tr, pp); // using interpolation scheme for temperpoint p
          m[] =tr.lambda*(tt - tr.tr_eq)/(L_h*d1);
        }
    }
  boundary({m});
}

// Simplified mass transfer inspired from Jie zhang's paper
// not use interpolation method, using averaging method
void zhang_model_1 (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  scalar f_v[];
  vertex scalar t_cor[];
  foreach_vertex()
    t_cor[] = (tr[] + tr[-1,0] + tr[0,-1] + tr[-1,-1])/4;
  boundary({t_cor});
  scalar t_av[];
  foreach()
    {
      m[] = 0;
      f_v[] = 1 - f[];
      coord n  = interface_normal( point, f) , s;
      double alpha  = plane_alpha (f[], n);
      plane_area_center(n, alpha, &s);
      double xc = 0, yc = 0;
      double delta_d = 1.5;// according to the paper, this value is between 1 and 2
      if(interfacial(point,f))
        {
          xc += s.x*Delta-n.x*delta_d*Delta;
          yc += s.y*Delta-n.y*delta_d*Delta;
          t_av[] = (t_cor[] + t_cor[1,0] + t_cor[0,1] + t_cor[1,1])/4;
          double d1 = fabs(n.x*xc + n.y*yc - alpha*Delta - Delta*0.5*(n.x+n.y))/sqrt(sq(n.x) + sq(n.y));// distance for temperature gradient
          m[] =tr.lambda*(t_av[] - tr.tr_eq)/(L_h*d1);
        }
    }
  boundary({m});
}

static inline double interp3 (Point point, coord p, scalar col) {
  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);
}

// Simplified mass transfer inspired from Jie zhang's paper
// using barycenter value in interfacial cell
void zhang_model_2 (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    {
      if(f[] == 1)
        tr[] = tr.tr_eq;
      f[] = clamp(f[], 0., 1.);
    }
  boundary({tr,f});

  scalar f_v[];
  vertex scalar t_cor[];
  foreach_vertex()
    t_cor[] = (tr[] + tr[-1,0] + tr[0,-1] + tr[-1,-1])/4;
  boundary({t_cor});
  face vector s[];
  foreach()
    {
      m[] = 0;
      f_v[] = 1 - f[];
      coord n  = facet_normal( point, f, s) , p, q, k;
      q.x = -n.x, q.y = -n.y;
      double alpha_new  = plane_alpha (f_v[], q);
      double alpha  = plane_alpha (f[], n); 
      line_center (q, alpha_new,f_v[],&p);// barycenter of vapor portion in interfacial cell
      line_center (n, alpha,f[],&k);// barycenter of liquid portion in interfacial cell
      double xc = 0, yc = 0;
      double xk = 0, yk = 0;
      if(f[]<1-1e-12&&f[]>1e-12)
        {
          xc += p.x*Delta;
          yc += p.y*Delta;
          xk += k.x*Delta;
          yk += k.y*Delta;
          coord pp;
          pp.x = xc, pp.y = yc;// mirror point in local corner neighboring cell
          double d1 = fabs(q.x*xc + q.y*yc - alpha_new*Delta - Delta*0.5*(q.x+q.y))/sqrt(sq(q.x) + sq(q.y));// distance for temperature gradient
          //double d1 = sqrt(sq(xk - xc) + sq(yk - yc));
          //double tt = interpolate_1 (point, tr, pp); // using interpolation scheme for temperpoint p
          double tt = interp3(point,p,tr);
          m[] = tr.lambda*(tt - tr.tr_eq)/(L_h*d1);// new method here
        }
    }
  boundary({m});
}

void Malan (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  double d[3][4];
  double eps[3][4];
  double gamma[3][4];
  double phi[3][4];
  double mm[3][4];
  //double m_sum[4];
  double tt[3][4]; 
  
  foreach()
  {
    coord p[3] = {{-1,1},{-1,0},{0,1}};
    coord q[3] = {{-1,-1},{-1,0},{0,-1}};
    coord r[3] = {{1,1},{1,0},{0,1}};
    coord s[3] = {{1,-1},{1,0},{0,-1}};
    coord n = interface_normal(point,f);
    normalize(&n);
    double alpha  = plane_alpha (f[], n) + 0.5*(n.x + n.y);  
    if(f[]>1e-12&&f[]<1-1e-12)
    {
      tt[0][0] = tr[-1,1];
      tt[1][0] = tr[-1,0];
      tt[2][0] = tr[0,1];
      tt[0][1] = tr[-1,-1];
      tt[1][1] = tr[-1,0];
      tt[2][1] = tr[0,-1];
      tt[0][2] = tr[1,1];
      tt[1][2] = tr[1,0];
      tt[2][2] = tr[0,1];
      tt[0][3] = tr[1,-1];
      tt[1][3] = tr[1,0];
      tt[2][3] = tr[0,-1];

      if(n.x<0&&n.y>0)
        {
          d[0][0] = fabs(-alpha + p[0].x*n.x + p[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][0] = fabs(p[0].x*n.x + p[0].y*n.y);
          gamma[0][0] = eps[0][0]/(sq(p[0].x) + sq(p[0].y));
          phi[0][0] = gamma[0][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[0][0] = phi[0][0]*(tt[0][0] - tr.tr_eq)/(Delta*d[0][0]);

          d[1][0] = fabs(-alpha + p[1].x*n.x + p[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][0] = fabs(p[1].x*n.x + p[1].y*n.y);
          gamma[1][0] = eps[1][0]/(sq(p[1].x) + sq(p[1].y));
          phi[1][0] = gamma[1][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[1][0] = phi[1][0]*(tt[1][0] - tr.tr_eq)/(Delta*d[1][0]);

          d[2][0] = fabs(-alpha + p[2].x*n.x + p[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][0] = fabs(p[2].x*n.x + p[2].y*n.y);
          gamma[2][0] = eps[2][0]/(sq(p[2].x) + sq(p[2].y));
          phi[2][0] = gamma[2][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
          mm[2][0] = phi[2][0]*(tt[2][0] - tr.tr_eq)/(Delta*d[2][0]);

          //m[] = mm[0][0] + mm[1][0] + mm[2][0];
        }
        
      if(n.x<=0&&n.y<=0)
        {
          d[0][1] = fabs(-alpha + q[0].x*n.x + q[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][1] = fabs(q[0].x*n.x + q[0].y*n.y);
          gamma[0][1] = eps[0][1]/(sq(q[0].x) + sq(q[0].y));
          phi[0][1] = gamma[0][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[0][1] = phi[0][1]*(tt[0][1] - tr.tr_eq)/(Delta*d[0][1]);
        
          d[1][1] = fabs(-alpha + q[1].x*n.x + q[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][1] = fabs(q[1].x*n.x + q[1].y*n.y);
          gamma[1][1] = eps[1][1]/(sq(q[1].x) + sq(q[1].y));
          phi[1][1] = gamma[1][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[1][1] = phi[1][1]*(tt[1][1] - tr.tr_eq)/(Delta*d[1][1]);
        
          d[2][1] = fabs(-alpha + q[2].x*n.x + q[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][1] = fabs(q[2].x*n.x + q[2].y*n.y);
          gamma[2][1] = eps[2][1]/(sq(q[2].x) + sq(q[2].y));
          phi[2][1] = gamma[2][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
          mm[2][1] = phi[2][1]*(tt[2][1] - tr.tr_eq)/(Delta*d[2][1]);

          m[] = mm[0][1] + mm[1][1] + mm[2][1];
        }
        

      if(n.x>=0&&n.y>=0)
        {
          d[0][2] = fabs(-alpha + r[0].x*n.x + r[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][2] = fabs(r[0].x*n.x + r[0].y*n.y);
          gamma[0][2] = eps[0][2]/(sq(r[0].x) + sq(r[0].y));
          phi[0][2] = gamma[0][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[0][2] = phi[0][2]*(tt[0][2] - tr.tr_eq)/(Delta*d[0][2]);
        
          d[1][2] = fabs(-alpha + r[1].x*n.x + r[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][2] = fabs(r[1].x*n.x + r[1].y*n.y);
          gamma[1][2] = eps[1][2]/(sq(r[1].x) + sq(r[1].y));
          phi[1][2] = gamma[1][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[1][2] = phi[1][2]*(tt[1][2] - tr.tr_eq)/(Delta*d[1][2]);
        
          d[2][2] = fabs(-alpha + r[2].x*n.x + r[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][2] = fabs(r[2].x*n.x + r[2].y*n.y);
          gamma[2][2] = eps[2][2]/(sq(r[2].x) + sq(r[2].y));
          phi[2][2] = gamma[2][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
          mm[2][2] = phi[2][2]*(tt[2][2] - tr.tr_eq)/(Delta*d[2][2]);

          m[] = mm[0][2] + mm[1][2] + mm[2][2];
        }

      if(n.x>0&&n.y<0)
        {
          d[0][3] = fabs(-alpha + s[0].x*n.x + s[0].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[0][3] = fabs(s[0].x*n.x + s[0].y*n.y);
          gamma[0][3] = eps[0][3]/(sq(s[0].x) + sq(s[0].y));
          phi[0][3] = gamma[0][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[0][3] = phi[0][3]*(tt[0][3] - tr.tr_eq)/(Delta*d[0][3]);
       
          d[1][3] = fabs(-alpha + s[1].x*n.x + s[1].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[1][3] = fabs(s[1].x*n.x + s[1].y*n.y);
          gamma[1][3] = eps[1][3]/(sq(s[1].x) + sq(s[1].y));
          phi[1][3] = gamma[1][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[1][3] = phi[1][3]*(tt[1][3] - tr.tr_eq)/(Delta*d[1][3]);
        
          d[2][3] = fabs(-alpha + s[2].x*n.x + s[2].y*n.y)/sqrt(sq(n.x) + sq(n.y));
          eps[2][3] = fabs(s[2].x*n.x + s[2].y*n.y);
          gamma[2][3] = eps[2][3]/(sq(s[2].x) + sq(s[2].y));
          phi[2][3] = gamma[2][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
          mm[2][3] = phi[2][3]*(tt[2][3] - tr.tr_eq)/(Delta*d[2][3]);

          m[] = mm[0][3] + mm[1][3] + mm[2][3];
        }
        
    }
    m[] *= tr.lambda/L_h;
  }
  boundary({m});
}

/*
void Malan (scalar tr, scalar f, scalar m, double L_h)
{
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary({f});

  double d[3][4];
  double eps[3][4];
  double gamma[3][4];
  double phi[3][4];
  double mm[3][4];
  double m_sum[4];


  coord p[3] = {{-1.5,1.5},{-1.5,0.5},{-0.5,1.5}};
  coord q[3] = {{-1.5,-1.5},{-1.5,-0.5},{-0.5,-1.5}};
  coord r[3] = {{1.5,1.5},{1.5,0.5},{0.5,1.5}};
  coord s[3] = {{1.5,-1.5},{1.5,-0.5},{0.5,-1.5}};
     
  foreach()
  {
    coord n = interface_normal(point,f);
    normalize(&n);

    double alpha  = plane_alpha (f[], n) + 0.5*(n.x + n.y); 
    d[0][0] = alpha + p[0].x*n.x + p[0].y*n.y;
    eps[0][0] = fabs(p[0].x*n.x + p[0].y*n.y);
    gamma[0][0] = eps[0][0]/(sq(p[0].x) + sq(p[0].y));
    phi[0][0] = gamma[0][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
    d[1][0] = alpha + p[1].x*n.x + p[1].y*n.y;
    eps[1][0] = fabs(p[1].x*n.x + p[1].y*n.y);
    gamma[1][0] = eps[1][0]/(sq(p[1].x) + sq(p[1].y));
    phi[1][0] = gamma[1][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
    d[2][0] = alpha + p[2].x*n.x + p[2].y*n.y;
    eps[2][0] = fabs(p[2].x*n.x + p[2].y*n.y);
    gamma[2][0] = eps[2][0]/(sq(p[2].x) + sq(p[2].y));
    phi[2][0] = gamma[2][0]/(gamma[0][0] + gamma[1][0] + gamma[2][0]);
    
    if(n.x<0&&n.y>0)
    {
      if(interfacial(neighborp(-1,-1),f)&&f[] == 0)
        mm[0][0] = phi[0][0]*(tr[] - tr.tr_eq)/(Delta*d[0][0]);
      else if(interfacial(neighborp(-1,0),f)&&f[] == 0)
        mm[1][0] = phi[1][0]*(tr[] - tr.tr_eq)/(Delta*d[1][0]);
      else if(interfacial(neighborp(0,1),f)&&f[] == 0)
        mm[2][0] = phi[2][0]*(tr[] - tr.tr_eq)/(Delta*d[2][0]);
    
    m[] = mm[0][0] + mm[1][0] + mm[2][0];
    }

    d[0][1] = alpha - q[0].x*n.x - q[0].y*n.y;
    eps[0][1] = fabs(q[0].x*n.x + q[0].y*n.y);
    gamma[0][1] = eps[0][1]/(sq(q[0].x) + sq(q[0].y));
    phi[0][1] = gamma[0][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
    d[1][1] = alpha - q[1].x*n.x - q[1].y*n.y;
    eps[1][1] = fabs(q[1].x*n.x + q[1].y*n.y);
    gamma[1][1] = eps[1][1]/(sq(q[1].x) + sq(q[1].y));
    phi[1][1] = gamma[1][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);
    d[2][1] = -alpha - q[2].x*n.x - q[2].y*n.y;
    eps[2][1] = fabs(q[2].x*n.x + q[2].y*n.y);
    gamma[2][1] = eps[2][1]/(sq(q[2].x) + sq(q[2].y));
    phi[2][1] = gamma[2][1]/(gamma[0][1] + gamma[1][1] + gamma[2][1]);     

    if(n.x<0&&n.y<0)
    {
      if(interfacial(neighborp(-1,-1),f)&&f[] == 0)
        mm[0][1] = phi[0][1]*(tr[] - tr.tr_eq)/(Delta*d[0][1]);
      else if(interfacial(neighborp(-1,0),f)&&f[] == 0)
        mm[1][1] = phi[1][1]*(tr[] - tr.tr_eq)/(Delta*d[1][1]);
      else if(interfacial(neighborp(0,-1),f)&&f[] == 0)
        mm[2][1] = phi[2][1]*(tr[] - tr.tr_eq)/(Delta*d[2][1]);
      m[] = mm[0][1] + mm[1][1] + mm[2][1];
    }


    d[0][2] = alpha - r[0].x*n.x - r[0].y*n.y;
    eps[0][2] = fabs(r[0].x*n.x + r[0].y*n.y);
    gamma[0][2] = eps[0][2]/(sq(r[0].x) + sq(r[0].y));
    phi[0][2] = gamma[0][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
    d[1][2] = alpha - r[1].x*n.x - r[1].y*n.y;
    eps[1][2] = fabs(r[1].x*n.x + r[1].y*n.y);
    gamma[1][2] = eps[1][2]/(sq(r[1].x) + sq(r[1].y));
    phi[1][2] = gamma[1][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);
    d[2][2] = alpha - r[2].x*n.x - r[2].y*n.y;
    eps[2][2] = fabs(r[2].x*n.x + r[2].y*n.y);
    gamma[2][2] = eps[2][2]/(sq(r[2].x) + sq(r[2].y));
    phi[2][2] = gamma[2][2]/(gamma[0][2] + gamma[1][2] + gamma[2][2]);

    if(n.x>0&&n.y<0)
    {
      if(interfacial(neighborp(1,1),f)&&f[] == 0)
        mm[0][2] = phi[0][2]*(tr[] - tr.tr_eq)/(Delta*d[0][2]);
      else if(interfacial(neighborp(1,0),f)&&f[] == 0)
        mm[1][2] = phi[1][2]*(tr[] - tr.tr_eq)/(Delta*d[1][2]);
      else if(interfacial(neighborp(0,1),f)&&f[] == 0)
        mm[2][2] = phi[2][2]*(tr[] - tr.tr_eq)/(Delta*d[2][2]);
      m[] = mm[0][2] + mm[1][2] + mm[2][2];
    }


    d[0][3] = alpha - s[0].x*n.x - s[0].y*n.y;
    eps[0][3] = fabs(s[0].x*n.x + s[0].y*n.y);
    gamma[0][3] = eps[0][3]/(sq(s[0].x) + sq(s[0].y));
    phi[0][3] = gamma[0][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
    d[1][3] = alpha - s[1].x*n.x - s[1].y*n.y;
    eps[1][3] = fabs(s[1].x*n.x + s[1].y*n.y);
    gamma[1][3] = eps[1][3]/(sq(s[1].x) + sq(s[1].y));
    phi[1][3] = gamma[1][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);
    d[2][3] = alpha - s[2].x*n.x - s[2].y*n.y;
    eps[2][3] = fabs(s[2].x*n.x + s[2].y*n.y);
    gamma[2][3] = eps[2][3]/(sq(s[2].x) + sq(s[2].y));
    phi[2][3] = gamma[2][3]/(gamma[0][3] + gamma[1][3] + gamma[2][3]);


    if(n.x>0&&n.y<0)
    {
      if(interfacial(neighborp(1,-1),f)&&f[] == 0)
        mm[0][3] = phi[0][3]*(tr[] - tr.tr_eq)/(Delta*d[0][3]);
      else if(interfacial(neighborp(1,0),f)&&f[] == 0)
        mm[1][3] = phi[1][3]*(tr[] - tr.tr_eq)/(Delta*d[1][3]);
      else if(interfacial(neighborp(0,-1),f)&&f[] == 0)
        mm[2][3] = phi[2][3]*(tr[] - tr.tr_eq)/(Delta*d[2][3]);
      m[] = mm[0][3] + mm[1][3] + mm[2][3];
    }

    m[] *= tr.lambda/L_h;
  }
  boundary({m});
}

*/


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
  face vector f_v[];
  foreach_face()
    {
      dd[] = 1 - f[];
      f_v.x[] = (dd[1] - dd[-1])/(2*Delta);
    }
  boundary((scalar *){f_v});

 /*
  foreach_face()
  {
    v_pc.x[] = (tr[1] - tr[-1])/(2*Delta);
  }
  boundary((scalar *){v_pc}); */
  foreach()
    {
      dd[] = 1 - f[];
      temp[] = 0;
      foreach_dimension()
        {
          if(f[]>1e-12&&f[]<1-1e-12)
            temp[] += 0.7*0.025*v_pc.x[]*f_v.x[]/L_h;
          else
            temp[] += 0;
      }
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

  scalar dirichlet_source_term[], volumic_metric[], b[];
  face vector diffusion_coef[];
  foreach()
    {
      volumic_metric[] = cm[];
      dirichlet_source_term[] = m[]*L_h/(tr.cp*tr.rho);
      b[] = 0;
    }
  boundary({volumic_metric,dirichlet_source_term,b});
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;
  boundary((scalar *){diffusion_coef});
  
  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term, beta = b, theta = volumic_metric);
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