// call reduced-gravity method
#define REDUCED 1
// adaptive method (changable method to 'adapt_limit.h')
#define ADAPT 1
// comparing with experimental results from chengsi
#define EXPERIMENT 1
// reference data from article of Araujo(2012)
#define PARTICLE 0

// axisymmetrical coordinate (to simulate 3D cases)
#include "axi.h" 

/********************* calculation cases ******************************/
#if EXPERIMENT
#include "experiment.h"
#endif

#if REFERENCE
#include "reference.h"
#endif
/***********************************************************************/

#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
// for large density ratio
#include "navier-stokes/conserving.h"
#include "tension.h"

#if REDUCED
  #include "reduced.h"
#endif
#include "adapt_limit.h"

// output
#include "view.h"
#include "vtk.h"

// record simulation statistics
#include "profiling.h"
#include "navier-stokes/perfs.h"

#if PARTICLE
#include "particle.h"
#include "scatter.h"
#endif

/**********************************************************/
#define Fr_0 1
#define _g 1 //in this case, consider U_0 = sqrt(gD_0)--->Fr_0^2 = U_0^2/(gD_0)

/* R_0 characteristic length (radias of lower tube) */
#define lambda_1  (2*(sqrt(1+2.44*pow(N_f,0.667))-1)/(2.44*pow(N_f,0.667))) 
// dimensionless liquid film thickness------refernce from chengsi's thesis tto guarantee the consistency with experiment
#define R_0 0.5
#define lambda  lambda_1*R_0 // liquid film
//#define r_0 (R_0 - lambda)// or (R_0 - lambda)
#define R_ex R_0*ratio
#define L_lower 8
#define L_upper 8
#define L0 L_lower+L_upper

/* grid infromation */
#define maxlevel 11
#define minlevel 8

// define two top walls
bid top_lower;
bid top_upper;

// no-slip boundary conditions on top walls
u.t[top_lower] = dirichlet(0);
u.t[top_upper] = dirichlet(0);

// right side for outflow
u.n[right] = neumann(0);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

/*there is no flow through the top and bottom boundary, 
otherwise the compatibility condition for the Poisson equation can be violated */
uf.n[bottom] = 0;
uf.n[top_lower] = 0;
uf.n[top_lower] = 0;

#define TIME 30

// geometry of Taylor-bubble
double sphere(double x, double y, coord center, double radius) {

  return ( sq(x ) + sq (y ) - sq (radius));

}

// Taylor bubble geometry
double geometry(double x, double y) {

  coord center;
  foreach_dimension()
  center.x = 0;
  double distance = 6;
  double s = sphere (x-(length_ratio+distance)*r_0, y, center, r_0);
  double left = x - distance*r_0;
  double right = -x + (length_ratio+distance)*r_0;
  double top = -y + r_0;
  double bottom=0;
  double c = -min(min(min(top,left),right),bottom);
  double geom = min(s, c);

  return geom;
}

int main() {
  size (L0);
  init_grid(1<<minlevel);
  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./N_f;
  mu2 = mu1/MUR;
  f.sigma = 1./Eo;
  TOLERANCE = 1e-3;

  #if REDUCED
    G.x -= _g;
  #endif
  run();  
}

void createFolder(const char* folder)

{
    char *str = NULL;
    str = (char *) malloc(sizeof(char) *1024);
    sprintf(str, "mkdir %s", folder); 
    system(str);   
}
#define dR_refine (2.*L0/(1 << minlevel))

event init (t = 0) {
// creat folders
createFolder("pressure");
createFolder("interface");
createFolder("tip");
createFolder("gfsview");
createFolder("velocity");
createFolder("field");
createFolder("shear");
createFolder("curvature");
createFolder("log_1");
createFolder("movie");
//createFolder("vortex");

// computational domain
  mask((x<=L_lower&&y>R_0)?top_lower:
  (x>L_upper&&y>(R_ex))?top_upper: none);
  
// mesh refine initallly for interested region
refine(x<=L_lower&&(y-R_ex<0)&&level<maxlevel);
refine(x>L_upper&&(y-R_ex<0)&&level<maxlevel);


fraction (f, geometry(x,y));
#if PARTICLE
double a[] = {8.2,7.5,9.5};
double b[] ={0.1,0.45,0.1};
// particle tracking
init_particles_random_circle(3, a,b, 0.001);
#endif
}

/**
We add the acceleration of gravity. */
#if !REDUCED
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= _g;
}
#endif

event logfile (t = 0; t <= TIME; t += 0.01) {

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 1024);
  sprintf (name, "log_1/out");
  static FILE * fp = fopen (name, "w");
  double xb = 0., yb = 0, vx = 0., vy = 0,  sb = 0., ke = 0;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:vx)  reduction(+:vy) reduction(+:sb) reduction(+:ke)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv;
    ke += dv()*rho(f[])*(sq(u.x[]) + sq(u.y[]));
  }
  
  fprintf (fp,"%g %g %g %g %g %g %g \n", t, sb, xb/sb, yb/sb, vx/sb, vy/sb, ke);
  fflush (fp);

}


#if ADAPT
event adapt (i++) {
  //adapt_wavelet ({f,u}, (double[]){1e-3,1e-2,1e-2}, maxlevel,minlevel);
  adapt_wavelet_leave_interface((scalar *){u},{f},(double[]){1e-1,1e-2,1e-2}, maxlevel, minlevel,0);
}
#endif


/**********************************************************/
#include "output_taylor.h"

// movie
event movies (t = 0; t <= TIME; t += 0.05) {
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  char legend_time[80];
  sprintf(legend_time, "t = %g s", t);
  
  clear();
// interface
  view (fov = 20, tx=-0.5,ty=0);
  //travelling (0, 0.5, fov = 20, tx = -0.5, ty = 0);
  travelling (0.5,TIME, tx=0,ty=0);
  draw_vof ("f",fc = {1,0,0}, lw = 2);
  squares("f",linear = true);
  #if PARTICLE
  scatter(loc = loc, s = 20, pc = {1, 0, 0});
  #endif
  draw_string(legend_time, 1, size = 30., lw = 2.);
  
  mirror({0,1}){
  draw_vof ("f",fc = {1,0,0}, lw = 2);
    }

  char name[100];
  sprintf (name, "movie/Taylor-bubble.mp4");
  save(name);
}