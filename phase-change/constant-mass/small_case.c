/*
note: to see the divergence velocity, not 'Devergence' in gfsview, I have defined several scalar filed.
*/
//#include "axi.h"
#include "centered-pc.h" // change the projection step
#include "two-phase-pc.h"
#include "tracer.h" // solving passive scalar advection 
#include "diffusion-pc.h"
#include "phase-change.h"
#include "tension.h"
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

//uf.n[left] = 0;
//uf.n[bottom] = 0;


int main() {
  size (1); // square domain 1m each side
  init_grid(1<<level); 

  rho1 = 1;
  rho2 = 0.001;
  //we need surface tension to keep droplet circle
  f.sigma = 0.0001;

  // poisson tolerance
  TOLERANCE = 1e-4;
  DT = 0.001;
  run();
}

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
  
}

/************************ define mass source *******************************/
scalar m[], div_pc[], delta_pc[];
scalar velocity[];     // velocity magnitude

event stability(i++)   // executed before vof and ns events
{
  mass_source(f,m);

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
    if(f[]>1e-12&&f[]<1-1e-12){
      coord n = interface_normal( point, f); // interface normal
      double alpha = plane_alpha(f[],n); // distance from original point to interface 
      alpha -= m[]*dt/(rho1*Delta)*sqrt(sq(n.x)+sq(n.y)); // phase-change shifting distance
      f[] = line_area(n.x,n.y,alpha); // cell volume fraction for n+1
      //double dc = 0.001*c[]; // (c[] - c.tr_eq)/dt
      //double df = 10*dc; // rho1*df*dv = m = dc/L_h -> df = dc/(L_h*dv)
      //f[] -= df;
      //c[] -= dc;
      //div_pc[] = df*((rho1)/(rho2)-1)/dt;
      div_pc[] = (f_old - f[])*((rho1)/(rho2)-1)/dt;
    }
    else
    {
      div_pc[] = 0;
    }
    
  }

  
//   delta_magnini( f, delta_pc );
//   foreach() 
//     div_pc[] = 0.05*delta_pc[]*(1./rho2+1./rho1);

//   smoother( div_pc, div_1 );

   
}

/*********************************** method to balance mass-weighted and volume-weighted velocity ********************************************************/
/*
This method is proposed by 'Direct numerical simulation of evaporating droplets', the gas expansion when droplet evaporation will cause difference between mass-weighted and volume-weighted velocity.
In order to make them balance, one new divergence-free source term replace the mass source in previous continuity equation. 
THis method was also used in 'Phase Change Heat Transfer Simulation for Boiling Bubbles Arising from a Vapor Film by the VOSET Method''
*/
/*********************************************************************************************/
face vector  s[], ull[], ugg[], uff[];
scalar div_new[]; 
face vector uu[];

event projection(i++){

  face_fraction (f, s); // calculating the face fraction of liquid by volume fraction 
  face vector n[];
  compute_normal(f,n);
  foreach_face()
  {
    s.x[] = 1 - s.x[];
    uu.x[] = (uf.x[]*s.x[]*rho1 + uf.x[]*(1-s.x[])*rho2)/(s.x[]*rho1+(1-s.x[])*rho2);
    ull.x[] = uu.x[] - div_pc[]*Delta*n.x[]*rho2/(s.x[]*rho1+(1-s.x[])*rho2); //Eq.38
    ugg.x[] = uu.x[] + div_pc[]*Delta*n.x[]*rho1/(s.x[]*rho1+(1-s.x[])*rho2)*(s.x[]/(1-s.x[]+1e-2)); //Eq.39
    uff.x[] = ull.x[]*s.x[]*rho1 + ugg.x[]*(1-s.x[])*rho2; //Eq.40
  }
  boundary((scalar *){ugg,ull,uff});

// calculating the divergence of new velocity to replace the Eq.19
  foreach()
  {
    div_new[] = 0;
    foreach_dimension()
      div_new[] += (uff.x[1] - uff.x[])/(cm[]*Delta);
  }
  boundary({div_new});

  foreach() 
    div_1[] = div_new[];
  boundary({div_1});

  smoother( div_1, div_new );
  smoother( div_new, div_1);
  smoother( div_1, div_new );

  project_pc(uf,p,alpha,dt,4,div_pc); // using new mass source term
  
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
      div_1[] += (u_l.x[1] - u_l.x[0])/Delta;
  }
  boundary({div_1}); 
  
  // solve pressure equation with phase change velocity divergence
  scalar divv[];
  foreach(){
    divv[] = div_pc[];
    divv[] /= dt;
  }
  poisson (p_new, divv, alpha, tolerance = TOLERANCE/sq(dt), nrelax = 2);
  
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


/*
event boundary_interface_domain(i++){
  foreach_boundary(right)
  {
    coord n = interface_normal(point,f);
    if(interfacial(point,f)&&(n.x>0))
      {
        p[] = -p_new[1];
        foreach_dimension()
          u.x[] = -u_l.x[1];
      }
  }

  foreach_boundary(top)
  {
    coord n = interface_normal(point,f);
    if(interfacial(point,f)&&(n.x>0))
      {
        p[] = -p_new[0,1];
        foreach_dimension()
          u.x[] = -u_l.x[0,1];
      }
  }


  foreach_boundary(left)
  {
    coord n = interface_normal(point,f);
    if(interfacial(point,f)&&(n.x<0))
      {
        p[] = -p_new[-1];
        foreach_dimension()
          u.x[] = -u_l.x[-1];
      }
  }

  foreach_boundary(bottom)
  {
    coord n = interface_normal(point,f);
    if(interfacial(point,f)&&(n.x<0))
      {
        p[] = -p_new[0,-1];
        foreach_dimension()
          u.x[] = -u_l.x[0,-1];
      }
  }
}
*/


