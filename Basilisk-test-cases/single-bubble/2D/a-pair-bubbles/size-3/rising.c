#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
//#include "tracer.h"
#include "tension.h"
#include "view.h"
#if 1
# include "lambda2.h"
#endif
#include "vtk.h"

#include "navier-stokes/perfs.h"
#if REDUCED
# include "reduced.h" //balance pressure with interfacial forces
#endif


double maxruntime = HUGE;
mgstats mgd;

#define Re 35
#define Eo 10
#define RHOR 1000.
#define MUR 100.

int LEVEL = 8;//grid level

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. The velocity on wall equals to 0->no slip boundary(Dirichlet),we don't know the real value but know the gradient of variable on the wall->slip boundary(Neummann)*/

u.t[bottom]  = dirichlet(0);

/**
We make sure there is no flow through the top and bottom boundary,
otherwise the compatibility condition for the Poisson equation can be
violated. */
//uf is the auxiliary face velocity,associated centered pressure pf

u.t[top] = dirichlet(0);
p[top] = dirichlet(0.);
pf[top] = dirichlet(0.);
uf.n[top] = neumann(0);

int main (int argc, char * argv[]) {
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    maxruntime = atoi (argv[2]);

  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  leverl=7~11. */
  size (4);

   rho1 = 1000.;
   rho2 = rho1/RHOR;
   mu1 = rho1*sqrt(0.98*0.5)*0.5/Re;
   mu2 = mu1/MUR;
   f.sigma = (rho1-rho2)*0.98*0.25/Eo;
 
    N = 1 << 7;
    init_grid(N);
    DT=0.01;
    TOLERANCE = 1e-4;
    run();
}


event init (t = 0) {
  /**
  The domain is a rectangle. We only simulate half the bubble. */
  
if (!restore (file = "dump")) {
  /**
  The bubble is centered on (0.5,0) and has a radius of 0.25. */
  fraction (f, min((sq(x-2.) + sq(y-0.5) +sq(z-2.) - sq(0.15)),
  ((sq(x-2.) + sq(y-1.1) +sq(z-2.) - sq(0.25)))));
}
}

/**
We add the acceleration of gravity. */

event acceleration (i++) {

  face vector av = a;
  foreach_face(y)
    av.y[] -= 0.98;
}


/**
A utility function to check the convergence of the multigrid
solvers. */

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0,
	    mg.nrelax);
}

/**
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

event logfile (t=0; t<=8.;t += 0.01) {
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
	  reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    zb += z*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    vbz += u.z[]*dv;
    sb += dv;
  }
fprintf (fout,
	   " %g %g %g %g %g %g %g %g "
	   "%g %d %d %d %ld %g %g\n", 
	   t, sb,
	   xb/sb, yb/sb, zb/sb,
	   vbx/sb, vby/sb, vbz/sb,
	   dt, mgp.i, mgpf.i, mgu.i,
	   grid->tn, perf.t, perf.speed);
  fflush (fout);
}



// Basilisk view 

event movies (t = 0; t <= 8.; t += 0.1) {

 char legend[1000];
  sprintf(legend, "t = %0.2g", t);
  clear();

// interface
  view (tx = -0.5,ty=-0.5);
  draw_vof ("f");
   draw_string(legend, 1, size = 30., lw = 2.);
    box (notics=true,lw=2);

save("interface.mp4");

}

// snapshot gerris(gfsview)
event snapshot (t = 0; t <= 8.1; t += 1)
{
 dump (file = "dump");
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump (name);
}

event snap (t = 0; t <= 8.1; t += 1)

{ scalar l2[], omegay[];
  #if 1
  lambda2 (u, l2);
  foreach()
    omegay[] = (u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
  #endif
  scalar omega[];
  vorticity(u,omega);
  char name[80];
  sprintf (name, "snapshot-%g.gfs", t);
  output_gfs (file = name, t = t, list = {f,u,p,l2,rho,mu});
}

event images(t=0; t+=0.5; t<=8.1){
scalar * list;
bool linear;
static FILE * fp = fopen ("paraview.vtk", "w");

output_vtk (list={f},N,fp,linear = true);

}

event adapt (i++) {
  adapt_wavelet ({f}, (double[]){5e-4}, LEVEL);
}

// adaptive mesh
#if MPI
event runtime (t += 0.01) {
  mpi_all_reduce (perf.t, MPI_DOUBLE, MPI_MAX);
  if (perf.t/60 >= maxruntime) {
    dump (file = "dump"); // so that we can restart
    return 1; // exit
}
}
#endif

