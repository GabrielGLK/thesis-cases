

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "view.h"

#define REDUCED 1
#if REDUCED
# include "reduced.h" 
#endif

#include "navier-stokes/perfs.h"
#include "maxruntime.h"

mgstats mgd;

#define Re 35
#define Eo 125
#define RHOR 1000.
#define MUR 100.

#define minlevel 6
#define maxlevel 11

#define MAXTIME 2.1
#define L0 2

p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
u.n[right] = neumann(0);

uf.n[bottom] = 0.;
uf.n[top] = 0.;

int main (){

  /**
  The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
  leverl=7~11. */
  size (L0);
  N = 1 << 7;
  init_grid(N);

  rho1 = 1000.;
  rho2 = 1./RHOR;
  mu1 = 10;
  mu2 = mu1/MUR;
  f.sigma = 1.96; 
  TOLERANCE = 1e-4;
  #if REDUCED
    G.x -= 0.98;
  #endif
    run();
}


event init (t = 0) {

  mask (y > L0/2 ? top : none);

  refine (sq(x-1) + sq(y ) - sq(0.75) < 0 && level < maxlevel);
  fraction (f, sq(x-1) + sq(y) - sq(0.5));

  CFL = 0.5;
}

/**
We add the acceleration of gravity. */
#if !REDUCED
event acceleration (i++) {

  face vector av = a;
  foreach_face(x)
    av.x[] -= 1;
}
#endif

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

event logfile (t=0; t<=MAXTIME;t += 0.01) {
// different benchmark quantities (position of centroid, circularity, mean rising velocity)
// note: the circularity should be calculated separately
  double xb = 0., vx = 0.,vy = 0., sb = 0.,yb = 0;
  foreach(reduction(+:xb) reduction(+:vx) reduction(+:sb) reduction(+:yb) reduction(+:vy)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv;
  }
  printf (" %g %g %g %g %g %g %g %g %g %g ", 
	  t, sb, -1., xb/sb, yb/sb,vx/sb, vy/sb,dt, perf.t, perf.speed);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  putchar ('\n');
  fflush (stdout);
}


// Basilisk view 

event movies (t = 0; t <= MAXTIME; t += 0.05) {

 char legend[1000];
  sprintf(legend, "t = %0.2g", t);
  clear();

// interface
  view (tx = -0.5,ty=0);
  draw_vof ("f",lc = {1,0,0},lw = 2);
  draw_string(legend, 1, size = 30., lw = 2.);
  
  mirror({0,1}){
    draw_vof ("f",lc = {1,0,0},lw = 2);
    }
 
save("interface.mp4");

}

event interface (t = 0; t <= MAXTIME; t += 1.) {

  char *outfile2 = NULL;
  outfile2 = (char *) malloc(sizeof(char) * 256);
  sprintf(outfile2, "interface-%g", t);
  FILE * fp_interface = fopen (outfile2, "w");
  output_facets (f, fp_interface);
  fclose(fp_interface);

}

// velocity vector field
event velocity(t+=1){
  char name[100];
 sprintf (name, "vprof-%g", t);
  FILE*fp = fopen (name, "w");
for (double y = 0; y <= L0/2; y += 0.05)
for (double x = 0; x <= L0; x += 0.05)
   fprintf (fp, " %g %g %g %g\n", x,y,
	     interpolate (u.x, x, y), interpolate (u.y, x, y));
fflush(fp);
}

event snap (t = 0; t <= MAXTIME; t += 1)

{ 
  scalar omega[];
  vorticity(u,omega);
  char name[80];
  sprintf (name, "snapshot-%g.gfs", t);
  output_gfs (file = name, t = t, list = {f,u,p,rho,mu});
}

// adaptive mesh
event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet ({f,u}, (double[]){1e-2,uemax,uemax,uemax}, maxlevel,minlevel);
}


