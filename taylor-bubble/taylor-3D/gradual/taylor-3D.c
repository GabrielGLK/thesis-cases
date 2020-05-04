/**********************************************************************/
//Taylor-bubble 3D test cases, compared with experimental data from chengsi

#define REDUCED 1
#define ADAPT 1
#define HEAT 0  // heat transfer
#define VERTICAL 1 //vertical tube
#define GRADUAL 1 //gradual tube
#define no_change 0 //vertical tube no-change
#define non_stick 0
#define MASK 1
#define popinet_trick  0
#define DEGREE 0
#define REFERENCE 0


#if REFERENCE
  #include "reference.h"
  #define RHOR 1000
  #define MUR 100
  #define length_ratio 10
  #define ratio 1
#endif 

#define TIME 25

#if VERTICAL
  #include "experiment.h"
#endif

/*******************************************************************************************/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

#if dimension == 3
# include "lambda2.h"
#endif

#if REDUCED
# include "reduced.h"
#endif

#include "view.h"
#include "interface_vtu_3D.h"
#include "navier-stokes/perfs.h"
/*******************************************************************************************/



#if HEAT
scalar T[];
scalar *tracers = {T};
face vector D[];
#endif

#define Fr_0 1
#define _g 1 //in this case, consider U_0 = sqrt(gD_0)--->Fr_0^2 = U_0^2/(gD_0)

/* R_0 characteristic length (radias of lower tube) */
#define lambda_1  (2*(sqrt(1+2.44*pow(N_f,0.667))-1)/(2.44*pow(N_f,0.667))) // dimensionless liquid film thickness------refernce from chengsi's thesis tto guarantee the consistency with experiment
#define R_0 0.5
#define R_ex R_0*ratio
#define L_lower 8
#define L_upper 8
#define L0 L_lower+L_upper
/* grid infromation */
#define maxlevel 10
#define minlevel 8
#define delta_ratio (ratio - 1)
#define delta_length delta_ratio*R_0

#if DEGREE
#define degree 45
#define theta degree*pi/180
#define delta_x delta_length/tan(theta)
#endif  

#if !DEGREE
#define delta_x 1
#define theta atan(delta_length/delta_x)  //gradual degree
#endif

// boundary conditions
u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);
u.n[right] = neumann(0);
pf[right] = dirichlet(0);


bid cylinder1;
u.t[cylinder1] = dirichlet(0.);
u.r[cylinder1] = dirichlet(0);
pf[cylinder1] = neumann(0.);
p[cylinder1] = neumann(0.);
uf.n[cylinder1] =dirichlet(0);



bid cylinder2;
u.t[cylinder2] = dirichlet(0.);
u.r[cylinder2] = dirichlet(0);
uf.n[cylinder2] = dirichlet(0);
pf[cylinder2] = neumann(0.);
p[cylinder2] = neumann(0.);

#if GRADUAL
bid gra;
u.t[gra] = dirichlet(0.);
u.r[gra] = dirichlet(0);
uf.n[gra] = dirichlet(0);
pf[gra] = neumann(0.);
p[gra] = neumann(0.);
#endif


int main() {
  size (L0);
  init_grid (1 << minlevel);
  rho1 = 1.;
  rho2 = 1./RHOR;
  mu1 = 1./N_f;
  mu2 = 1./(MUR*N_f);
  f.sigma = 1./Eo;
  TOLERANCE = 1e-3;
#if REDUCED
  G.x -= 1;
  Z.x = 1;
#endif
  run();
}


double sphere(double x, double y, double z, coord center, double radius) {
  return ( sq(x - center.x) + sq (y - center.y) + sq (z - center.z)
	   - sq (radius));
}


double geometry(double x, double y, double z) {

  coord center;
  foreach_dimension()
    center.x = 0;
  double distance = 18;
  double s = sphere (x-(length_ratio+distance)*r_0, y-L0/2, z-L0/2, center, r_0);
  double left = x - distance*r_0;
  double right = -x + (length_ratio+distance)*r_0;
  double cylind = sq(y-L0/2)  +sq(z-L0/2)- sq(r_0);
  double c = max(-min(left,right),cylind);

  double geom = min(s, c);

  return geom;
}


void createFolder(const char* folder)
{
    char str[128];
    sprintf(str, "mkdir %s", folder); 
    system(str);   
}


event init (t = 0) {

createFolder("interface");

#if no_change
mask(sq(y-L0/2) + sq(z-L0/2) > sq(R_0)? cylinder1::none);
#endif



#if GRADUAL
  mask(x<=L_lower&&sq(y-L0/2) + sq(z-L0/2) > sq(R_0)? cylinder1: 
  (L_lower<x&&x<=(L_lower+delta_x)&&sq(y-L0/2) + sq(z-L0/2) - sq(R_0)>(x-L_lower)*tan(theta))?gra:
  (x>(L_lower+delta_length)&&sq(y-L0/2) + sq(z-L0/2) > sq(R_ex))? cylinder2:
  none);
 #endif

  refine(8*r_0<x&&x<(length_ratio+5)*r_0 &&level<maxlevel);
  fraction (f, geometry(x,y,z));

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
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}

/*
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

event logfile (t=0; t<=TIME;t += 0.01) {
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



#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f}, (double[]){1e-2}, maxlevel,minlevel);
}
#endif

/*************************************output******************************************************/


// output video file
event movies (t = 0; t <= TIME; t += 0.2) {
  
  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

  char legend[80];
  sprintf(legend, "t = %0.2g", t);
  clear();

// interface
  view (tx = -0.5,ty=-0.5);
  draw_vof ("f");
  //foreach_face()
    //u.x[] += y;
  //for (double iv = 0; iv <= 0.5; iv += .1) 
    //iso_contour(u.x, iv);
  draw_string(legend, 1, size = 30., lw = 2.);

  save("interface.mp4");
  }


// this is file transfer using in Paraview, opacity
event save_vtu_surface (t = 0; t <= TIME; t += 0.1){
foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

char *name = NULL;
name = (char *) malloc(sizeof(char) * 256);
FILE *fp;
sprintf (name, "interface/paraview--%d-%.6f.vtu",pid(),t);
fp = fopen (name, "w");
output_vtu (f,fp);
fclose(fp);

}


event save_vtu_surface (t = 0; t <= TIME; t += 0.1){
foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f});

char *name = NULL;
name = (char *) malloc(sizeof(char) * 256);
FILE *fp;
sprintf (name, "interface/fied--%d-%.6f.vtu",pid(),t);
fp = fopen (name, "w");
output_vtu_ascii_foreach ({f,p,rho,u.x,u.y},{mu},N,fp,true);
fclose(fp);

}


