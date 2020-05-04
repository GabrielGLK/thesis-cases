//#include "axi.h"
#define JACOBI 1
#include "centered-pc.h" // change the projection step
#include "two-phase-pc.h"
#include "phase-change.h"
#include "tension.h"
#include "navier-stokes/conserving.h" // This file implements momentum-conserving VOF advection of the velocity components for the two-phase Navier-Stokes solver.
#include "navier-stokes/perfs.h" // record velocity/pressure statistic 
#include "view.h"

#define REDUCED 1
#if REDUCED
# include "reduced.h"
#endif

#define level 8

double maxruntime = HUGE;

#define L0 0.0788

#define T_sat 373.15// could not be zero, because in Lee model for denominator
#define T_sup 5 // difference between wall and saturated temperatures
#define T_wall (T_sat + T_sup)

#define cp_1 400
#define cp_2 200
#define lambda_1 40
#define lambda_2 1
#define D_L lambda_1/(rho1*cp_1)
#define D_V lambda_2/(rho2*cp_2)
#define L_h 1e4

face vector u_l[];    // liquid face velocity
face vector uf_new[];
scalar div_1[];       // one-fluid velocity divergence
scalar div_2[];       // liquid velocity divergence
scalar p_new[];

/*********************************** temperature tracers **************************************/
scalar T_V[], T_L[], *tracers = {T_V, T_L}; // vapor and liquid temp and use tracer.h to advect them

/************************************ boundary conditions **************************************/
T_V[bottom] = dirichlet(T_wall);

p[top] = dirichlet(0.);
u.n[top] = neumann(0.);
p_new[bottom] = dirichlet(0.);

uf.n[left] = 0;
uf.n[right] = 0;
uf.n[bottom] = 0;
mgstats mgd;


int main() {
  size (3*L0); // square domain 1m each side
  init_grid(1<<level); 
  rho1 = 200;
  rho2 = 5;
  mu1 = 0.1;
  mu2 = 0.005;
  f.sigma = 0.1;
  TOLERANCE = 1e-4;

  #if REDUCED
    G.y = -9.81;
  #endif
  run();
}

#define H0 L0/128
#define plane(x, y, H) (y-H0*cos(2*pi*(x-L0)/(L0)) - 4*H0)
#define char_L sqrt(f.sigma/((rho1-rho2)*9.81))
#define rhoo (1/rho2 - 1/rho1)

event init (t = 0) {
  mask(x>L0/2?right:none);
  fraction (f, plane(x, y, H0));
 
  foreach(){
    T_V[] = (1-f[])*(T_wall-T_sup/(H0*cos(2*pi*(x-L0)/(H0*128)) + 4*H0)*y);
    T_L[] = f[]*T_sat;
  }
  foreach_face(){
    uf.x[] = 0.;
}
  boundary({T_L, T_V, uf});

  //CFL = 0.2;
  }

/************************ define mass source *******************************/
scalar velocity[];     // velocity magnitude
scalar m_dot[];
scalar div_pc[];

event stability(i++)   // executed before vof and ns events
{
  T_L.rho = rho1;
  T_V.rho = rho2;
  T_L.lambda = lambda_1;
  T_V.lambda = lambda_2;
  for(scalar t in tracers)
    {
      t.tr_eq = T_sat;
      lee_mass(f,t,m_dot,L_h);
      mass_diffusion(div_pc,f,m_dot); //Hardt method
    }
  

  // compute velocity magnitude
  foreach()
    velocity[] = sqrt(sq(u.x[]) + sq(u.y[]));
  boundary({velocity});
}


event vof(i++)
{

/******************* step-2: shift interface accounting for phase-change *******************/  
  foreach()
  {
    double f_old = f[];
   /* Note that VOF function is clipped when the interface displacement extends beyond the cell boundary.
   * This could be solved by addig the clipped value to neighbouring cells. */
if(interfacial(point,f))
      {coord n = interface_normal( point, f); // interface normal
      double alpha = plane_alpha(f[],n); // distance from original point to interface 
      alpha -= m_dot[]*dt/(rho1*Delta)*sqrt(sq(n.x)+sq(n.y)); // phase-change shifting distance
      f[] = line_area(n.x,n.y,alpha); // cell volume fraction for n+1
  }
  }
  boundary({f});
}

event tracer_diffusion(i++){
  T_V.tr_eq = T_sat;
  T_L.tr_eq = T_sat;
  T_L.D = D_L;
  T_V.D = D_V;
  T_L.rho = rho1;
  T_V.rho = rho2;
  T_L.cp = cp_1;
  T_V.cp = cp_2;
  T_L.inverse = false;
  T_V.inverse = true;
  for(scalar t in tracers)
    heat_source (t, f, m_dot, L_h);
}

// first approximation projection method after prediction of face velocity to compute advection velocity, this step also considers mass source due to projection step
event advection_term (i++,last)
{
  if (!stokes) {
    mgpf = project_pc (uf, pf, alpha, dt/2., mgpf.nrelax,div_pc,rhoo);
  }
}


event projection(i++){
  
  project_pc(uf,p,alpha,dt,4,div_pc,rhoo); // projection step for entire domain
  
  centered_gradient (p, g);

  correction (dt);

  
  // copie velocity to u_l
  foreach_face()
    u_l.x[] = uf.x[];
  boundary((scalar *){u_l});

  // compute de divergence of u_l (should be exactly equal to div_pc)
  foreach(){
    div_1[] = 0;
    foreach_dimension()
      div_1[] += (u_l.x[1] - u_l.x[0])/(cm[]*Delta);
  }
  boundary({div_1}); 
  
  // solve pressure equation with phase change velocity divergence
  scalar divv[];
  foreach(){
    divv[] = div_pc[]*rhoo;
    divv[] /= dt;
  }
  poisson (p_new, divv, alpha, tolerance = TOLERANCE/sq(dt), nrelax = 4);
  
  // add phase change correction to u_l
  foreach_face()
    {
      uf_new.x[] = -dt*alpha.x[]*face_gradient_x (p_new, 0);
      u_l.x[] = uf.x[] + uf_new.x[];
    }
  boundary ((scalar *){u_l, uf_new});

  // compute again velocity divergence of u_l (should be zero)
  foreach(){
    div_2[] = 0;
    foreach_dimension()
      div_2[] += (u_l.x[1] - u_l.x[0])/Delta;
  }
  boundary({div_2});   

}

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0,
	    mg.nrelax);
}


event logfile (t = 0; t <= 2.1; t += 0.01) {

  double xb = 0., vx = 0.,vy = 0., sb = 0.,yb = 0., nu = 0.;
  foreach(reduction(+:xb) reduction(+:vx) reduction(+:sb) reduction(+:yb) reduction(+:vy) reduction(+:nu)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    nu += char_L/T_sup*(T_wall-(1-f[])*(T_wall-T_sup/(4*H0)*(y)))/y*dv();
    sb += dv; 
  }
  printf (" %g %g %g %g %g %g %g %g %g %g %g ", 
	  t, sb, -1., 2*nu/L0, xb/sb, yb/sb,vx/sb, vy/sb,dt, perf.t, perf.speed);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  putchar ('\n');
  fflush (stdout);
}

// output


event snap (t=0;t += 0.01)
 {
   char name[80];
   sprintf (name, "snapshot-%g.gfs",t);
   output_gfs (file = name);
 }


event movies (t = 0; t <= 2.1; t += 0.01) {
scalar T[];
  foreach()
    T[] = f[]*T_L[] + (1. - f[])*T_V[];
  boundary ({T});
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  char legend[1000];
  sprintf(legend, "t = %0.2g", t);
  clear();
// interface
  view (tx = 0,ty=-0.5);
  draw_vof ("f",lc = {0,0,0},lw = 2);
   draw_string(legend, 1, size = 30., lw = 2.);
    box (notics=true,lw=2);
     mirror({-1}){
    draw_vof("f",lc = {0,0,0},lw = 2);
    draw_string(legend, 1, size = 30., lw = 2.);
     box (notics=true,lw=2);
}
save("interface.mp4");
  }



