#define JACOBI 1 // relaxation method in Poisson equation
#include "../phase_change_code/embed-pc.h"
#include "../phase_change_code/centered-pc.h" // change the projection step
#include "../phase_change_code/two-phase-pc.h"
#include "../phase_change_code/phase-change.h"
//#include "tension-pc.h"
#include "navier-stokes/conserving.h" // This file implements momentum-conserving VOF advection of the velocity components for the two-phase Navier-Stokes solver.
#include "navier-stokes/perfs.h" // record velocity/pressure statistic 
#include "view.h"

// grid level
#define level 7

double maxruntime = HUGE;

face vector u_l[];    // liquid face velocity
face vector uf_new[];
scalar div_1[];       // liquid velocity divergence
scalar div_2[];       // liquid velocity divergence
scalar p_new[];
scalar c[], * tracers = {c};           // normalized temperature

mgstats mgd;// using diffusion solver

/************************************ boundary conditions **************************************/
// outflow boundary conditions
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
p_new[top] = dirichlet(0.);
u.n[top] = neumann(0.);
u_l.n[top] = neumann(0.);

p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
p_new[right] = dirichlet(0.);
u.n[right] = neumann(0.);
u_l.n[right] = neumann(0.);

uf.n[left] = 0;
uf.n[bottom] = 0;


int main() {
  size (1); // square domain 1m each side
  init_grid(1<<level); 

  rho1 = 1;
  rho2 = 0.01;
  //we need surface tension to keep droplet circle
  //f.sigma = 0.;

  // poisson tolerance
  TOLERANCE = 1e-4;
  DT = 0.001;
  run();
}


#define rhoo (1/rho2 - 1/rho1)
event init (t = 0) {  

  // geometry - static droplet with 0.23m radius
  fraction (f,  sq(0.23) - sq((x)) - sq((y)));
  
  // temperature equals 1 in liquid and 0 in gas
//   foreach()
//     c[] = f[];
//   boundary({c});
  
  // making superheated zone slightly smaller
  // for a soft start of phase change
  fraction (c,  sq(0.2) - sq((x)) - sq((y)));

  CFL = 0.4;
  
}

/************************ define mass source *******************************/
scalar m[], div_pc[], delta_pc[];
scalar velocity[];     // velocity magnitude

event stability(i++)   // executed before vof and ns events
{
  mass_source(f,m);
  mass_diffusion(div_pc,f,m); //Hardt method
  // compute velocity magnitude
  foreach()
    velocity[] = sqrt(sq(u.x[]) + sq(u.y[]));
  boundary({velocity});

}
scalar temp1[], temp2[];
event vof(i++)
{
/******************* step-2: shift interface accounting for phase-change *******************/  
  foreach()
  {
    double f_old = f[];
   /* Note that VOF function is clipped when the interface displacement extends beyond the cell boundary.
   * This could be solved by addig the clipped value to neighbouring cells. */
    if(f[]>1e-12&&f[]<1-1e-12){
      coord n = interface_normal( point, f); // interface normal
      double alpha = plane_alpha(f[],n); // distance from original point to interface 
      alpha -= m[]*dt/(rho1*Delta)*sqrt(sq(n.x)+sq(n.y)); // phase-change shifting distance
      f[] = line_area(n.x,n.y,alpha); // cell volume fraction for n+1
    }
}
  boundary({f});
}


event tracer_diffusion (i++)
{
  const face vector D[] = {0.001, 0.001};
  diffusion (c, dt, D );
}


event projection(i++){

  project_pc(uf,p,alpha,dt,4,div_pc,rhoo); 
  centered_gradient (p, g);
  correction (dt);
  
  // solve pressure equation with phase change velocity divergence
  scalar divv[];
  foreach(){
    divv[] = div_pc[]*rhoo;
    divv[] /= dt;
  }
  poisson (p_new, divv, alpha, tolerance = TOLERANCE/sq(dt), nrelax = 2);
  
  // copie velocity to u_l
  foreach_face()
    u_l.x[] = uf.x[];
  boundary((scalar *){u_l});

  // compute de divergence of u_l (should be exactly equal to div_pc)
  foreach(){
    div_1[] = 0;
    foreach_dimension()
      div_1[] += (u_l.x[1] - u_l.x[0])/(Delta*cm[]);
  }
  boundary({div_1}); 

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
      div_2[] += (u_l.x[1] - u_l.x[0])/(Delta*cm[]);
  }
  boundary({div_2});  

}


// event logfile (i=0;i<=10;i += 1) {
event logfile (t=0;t<1;t += 0.1) {
  double xb = 0., vb = 0., sb = 0., sm = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb) reduction(+:sm)) { // integration
    double dv = f[]*dv(); // the droplet volume in each cell
    sb += dv; // sum volume for all grids including liquid
    double dm = delta_pc[]*dv(); 
    sm += dm;
  }
  double ve = pi*sq(0.23-0.05*t); 
  printf ("%g %g %g %g %g",
          t, sb*4, sm*4, ve, dt ); // output real droplet volume

  putchar ('\n');
  fflush (stdout);
}

// event snap (i=0;i<=10;i += 1)
// {
//   char name[80];
//   sprintf (name, "snapshot-%i.gfs",i);
//   output_gfs (file = name);
// }
event snap (t=0;t += 0.1)
{
  char name[80];
  sprintf (name, "snapshot-%g.gfs",t);
  output_gfs (file = name);
}


event movies (t = 0; t += 0.01) {

  char legend[100];
  sprintf(legend, "t = %g", t);
  clear();
// interface
  view (tx = -0.5,ty=-0.5);
  draw_vof ("f",lc = {0,0,0},lw = 2);
  squares("m");
  draw_string(legend, 1, size = 30., lw = 2.);
  box (notics=true,lw=2);
     
  save("interface.mp4");
  }

