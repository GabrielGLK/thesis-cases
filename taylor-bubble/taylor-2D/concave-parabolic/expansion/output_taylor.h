/*
The output file for 2D cylinderical coordinate Taylor bubble rising, as follows:
1). bubble interface and cell evolution with time
2). the distance between bubble nose and expansion after bubble through expansion
3). coordinate of bubble nose and tail evolution with time
4). bubble frontal radius 
5). velocity field 
6). wall shear stress
7). pressure ---- middle liquid film (for bubble body)
8). gfsview file
9). output_vtk (for other scalar or vector field)
*/


/******************************** output **********************************************/
#include "fracface.h"


event interface (t += 0.1) {

// interface
    char *outfile1 = NULL;
    outfile1 = (char *) malloc(sizeof(char) * 256);
    sprintf(outfile1, "interface/interface-%.1f",t);
    FILE * fp1 = fopen (outfile1, "w");
    output_facets (f, fp1);
    fclose(fp1);
// cell 

    char *outfile = NULL;
    outfile = (char *) malloc(sizeof(char) * 1024);
    sprintf(outfile, "interface/cell-%.1f",t);  
    FILE * fp = fopen (outfile, "w");
    output_cells (fp);
    fclose (fp);

/*
    vertex scalar phi[];
    foreach_vertex()
        phi[] = -geometry(x,y);
    scalar c[];
    face vector sf[];
    fractions (phi,c);
    face_fraction (c,sf);

    scalar pos[];
    position (c, pos, {1,0});
    double max = statsf(pos).max;
    double min = statsf(pos).min;

    scalar pos_1[];
    position (c, pos_1, {0,1});
    double max_1 = statsf(pos_1).max;
    double min_1 = statsf(pos_1).min;

    char *name= NULL;
    name = (char *) malloc(sizeof(char) * 1024);
    
    
    for (x = min; x <= max; x += 0.1 )
        for (y = min_1; y <= max_1; y += 0.1)
            {
                sprintf (name, "interface/face-%g", t);
                static FILE * fp = fopen (name, "w");
                foreach_face()
                    fprintf(fp,"%g %g %g \n", x, y, sf.x[]);
                fclose (fp);
                }
                */
scalar neg[], pos[];

    foreach() {
        pos[] = f[];
        neg[] = (1. - f[]);
  }

char *name = NULL;
name = (char *) malloc(sizeof(char) * 1024);
sprintf(name, "interface/face-%.1f",t);  
FILE * fp2;
fp2 = fopen (name, "w");
foreach()
    fprintf(fp2,"%g %g %4.2f %4.2f %4.2f \n", x, y, pos[], neg[], f[]);
fclose (fp2);

}


// tip velocity
/*
event tip (t += 0.01) {
    scalar pos[];
    position (f, pos, {1,0});
    double max = statsf(pos).max;
    double min = statsf(pos).min;
    char name[80];
    sprintf (name, "tip/tip");
    static FILE * fp = fopen (name, "w");
    fprintf (fp, "%g %g %g\n", t, max,min);
    
    fflush (fp);
}
*/

// liquid film
event liquid (t += 0.01) {
    scalar pos[];
    position (f, pos, {1,0},{8,0});
    double max = statsf(pos).max;
    double min = statsf(pos).min;
    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "tip/liquid");
    static FILE * fp = fopen (name, "w");
    fprintf (fp, "%g %g %g\n", t, max,min);
    
    fflush (fp);
}

// frontal radius
event front_radius(t += 0.1){
    vector h[];
    heights (f, h); // return height 
    
    scalar pos[];
    position (f, pos, {1,0});
    double max = statsf(pos).max;
    scalar pos1[];
    position (f, pos1, {0,1});
    double max_1 = statsf(pos1).max;

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "curvature/curvature");
    static FILE * fp = fopen (name, "w");
    Point point = locate (max, 0);
    double a = height_curvature(point, f, h); // return kappa
    fprintf (fp, "%g %g %g %g %g\n", t, 1/fabs(a),max,max_1, 2/(fabs(a)*max_1));
    fflush (fp);
}

// position tail center

event particle1(t += 0.1){

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "curvature/point-1");
    static FILE * fp = fopen (name, "w");
    fprintf (fp, "%g %g %g %g %g %g\n", t,loc[0].x,loc[0].y, 
    interpolate (u.x, loc[0].x, loc[0].y), interpolate (u.y, loc[0].x, loc[0].y), interpolate (p, loc[0].x, loc[0].y));
    fflush (fp);
}

event particle2(t += 0.1){

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "curvature/point-2");
    static FILE * fp = fopen (name, "w");
    fprintf (fp, "%g %g %g %g %g %g\n", t,loc[1].x,loc[1].y, 
    interpolate (u.x, loc[1].x, loc[1].y), interpolate (u.y, loc[1].x, loc[1].y), interpolate (p, loc[1].x, loc[1].y));
    fflush (fp);
}

event particle3(t += 0.1){

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "curvature/point-3");
    static FILE * fp = fopen (name, "w");
    fprintf (fp, "%g %g %g %g %g %g\n", t,loc[2].x,loc[2].y, 
    interpolate (u.x, loc[2].x, loc[2].y), interpolate (u.y, loc[2].x, loc[2].y), interpolate (p, loc[2].x, loc[2].y));
    fflush (fp);
}


// position tail center
/*
event head_center(t += 0.1){

    char name[80];
    sprintf (name, "curvature/point_head");
    static FILE * fp = fopen (name, "w");
    fprintf (fp, "%g %g %g %g %g\n", t,loc[1].x,loc[1].y, 
    interpolate (u.x, loc[1].x, loc[1].y), interpolate (u.y, loc[1].x, loc[1].y) );
    fflush (fp);
}
*/


/*
// vortex
event vortex_1 (t += 0.1){
    scalar pos[];
    position (f, pos, {1,0});
    double max = statsf(pos).max;
    double min = statsf(pos).min;
    scalar vortex[];
    position(f, vortex, {1,0},{min,0});
    char name[80];
    sprintf (name, "vortex/vortex-%.1f",t);
    static FILE * fp = fopen (name, "w");
    fprintf (fp, "%g %g \n", x, statsf(vortex).max);
    fclose(fp);
}
*/

// output velocity distribution
event velocity (t+=0.1){

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "velocity/vprof-%.1f",t);
    FILE*fp1 = fopen (name, "w");


    for (double y = 0; y < R_ex; y += 0.02)
    for (double x = 0; x < L0; x += 0.1)
        fprintf (fp1, " %g %g %g %g %g\n", t, x,y,
            interpolate (u.x, x, y), interpolate (u.y, x, y)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
        fclose(fp1);
  }


event velocity_center (t+=0.1){

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "velocity/vcenter-%.1f",t);
    FILE*fp1 = fopen (name, "w");

    double vx = 0., vy = 0,  sb = 0.;
  foreach(reduction(+:vx)  reduction(+:vy) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    vx += u.x[]*dv;
    vy += u.y[]*dv;
    sb += dv;
  }

vector uu[];
foreach_face(){
  uu.x[] = u.x[] - vx/sb*fm.x[];
  uu.y[] = u.y[] - vy/sb*fm.y[];
}
boundary((scalar *){uu});


    for (double y = 0; y < R_0; y += 0.02)
    for (double x = 0; x < L0; x += 0.1)
        fprintf (fp1, " %g %g %g %g %g\n", t, x,y,
            interpolate (uu.x, x, y), interpolate (uu.y, x, y)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
        fclose(fp1);
  }


/*
event parameter (t+=0.1){

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "curvature/parameter-%.1f",t);
    FILE*fp1 = fopen (name, "a");

    scalar omega[];
    foreach()
        omega[] = (u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/(2.*Delta);
    boundary ({omega});

    //for (double y = 0; y <= 0.2; y += 0.02)
    for (double x = 0; x <= L0; x += 0.1)
        fprintf (fp1, "%g %g %g\n", t,x,interpolate (omega, x, 0.1)); 
    fflush(fp1);
  }
*/

// gfsview
event gfsview (t = 0; t <= TIME; t += 0.5)
{
    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    FILE *fp;
    sprintf (name, "gfsview/snapshot-%.1f.gfs",t);
    fp = fopen (name, "w");
    output_gfs (fp);
}

/*
event gfsview (t = 0; t <= TIME; t += 0.1)
{
  char name[80];
  sprintf(name,"gfsview2D -s 1.gfv | gfsview/snap-%.1f",t);
  static FILE * fp = popen (name, "w");
  output_gfs (fp, t = t);
  fprintf (fp, "Save { format = Gnuplot}");
}
*/

// output vtk with field

event save_vtu_surface (t = 0; t <=  TIME; t += 0.5){

    foreach()
        f[] = clamp(f[], 0., 1.);
    boundary ({f});

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    FILE *fp;
    sprintf (name, "field/paraview-%.1f.vtk",t);
    fp = fopen (name, "w");
    output_vtk ({f,p,rho,mu,u.x,u.y},N, fp, true);
    fflush(fp);
}


event shear (t+=0.1){

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "shear/shear-%.1f",t);
  FILE*fp1 = fopen (name, "w");
  scalar pos_1[];
  position (f, pos_1, {1,0});
  double max_1 = statsf(pos_1).max;
  double min_1 = statsf(pos_1).min;
  scalar shear[];
  foreach()
    shear[] = mu.y[]*(u.x[0,1] - u.x[0,-1])/(2*Delta);
  boundary ({shear});

  for (double y = (R_ex-0.002); y < R_ex; y += 0.002)
  for (double x = min_1-0.2; x <= max_1+0.2; x += 0.1)
    fprintf (fp1, " %g %g %g %g %g\n", x,y,
	    interpolate (u.x, x, y), interpolate (u.y, x, y),
      interpolate (shear, x, y)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
    fclose(fp1);
  }
  



  event shear2 (t+=0.1){

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "shear/shear1-%.1f",t);
  FILE*fp1 = fopen (name, "w");
  scalar pos_1[];
  position (f, pos_1, {1,0});
  double max_1 = statsf(pos_1).max;
  double min_1 = statsf(pos_1).min;
  scalar shear[];
  foreach()
    shear[] = mu.y[]*(u.x[0,1] - u.x[0,-1])/(2*Delta);
  boundary ({shear});

  for (double y = 0.498; y < R_0; y += 0.002)
  for (double x = min_1-0.2; x <= max_1+0.2; x += 0.1)
    fprintf (fp1, " %g %g %g %g %g\n", x,y,
	    interpolate (u.x, x, y), interpolate (u.y, x, y),
      interpolate (shear, x, y)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
    fclose(fp1);
  }



// pressure --- middle line of liquid film

  event pressure (t+=0.1){

  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;

  char *name = NULL;
  name = (char *) malloc(sizeof(char) * 256);
  sprintf (name, "pressure/pressure-%.1f", t);
  FILE*fp1 = fopen (name, "w");
  double b = max; 
  double a = 0;
  for (double x = 0; x <= L0; x += 0.1)
  {
  if (x<8)
    {a = (R_0 - b)/2 + r_0;}
  else
    {a = (R_ex -b)/2 + r_0;}
    fprintf (fp1, " %g %g %g\n", x,(a-r_0)*2,
	    interpolate (p, x, a)); // generate velocity vector field, at least 12 grids space in liquid thin film(\delta)----reference from "Co-current flow effects on a rising Taylor bubble"
    
  }
  fclose(fp1);
  }

event pressure (t+=0.1){

    scalar pos[];
    position (f, pos, {0,1});
    double max = statsf(pos).max;

    scalar pos_1[];
    position (f, pos_1, {1,0});
    double max_1 = statsf(pos_1).max;
    double min_1 = statsf(pos_1).min;

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "pressure/pressure-%.1f",t);
    FILE*fp = fopen (name, "w");
    
    for (double y = (max+(R_0 - max)/2); y < (max+(R_0 - max)/2) + 0.0001; y += 0.0001)
    for (double x = min_1; x <= max_1; x += 0.1)
        fprintf (fp, " %g %g %g\n", x,y,interpolate (p, x, y)); 
    fclose(fp);
  }



event mid_expansion(t+=0.1){
    scalar shear[];
    foreach()
        shear[] = mu.y[]*(u.x[0,1] - u.x[0,-1])/(Delta);
    boundary ({shear});
    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 256);
    sprintf (name, "pressure/mid-prop-%g-%g-%g-%g-%g",l_taylor,RHOR,MUR,N_f,Eo);
    FILE*fp = fopen (name, "a");
    fprintf (fp, " %g %g %g %g\n", t,interpolate (p, 8, 0), interpolate(u.x,8,0),interpolate(shear,8,0));
    fflush(fp);


}


/*
event waveprobe (t+=0.1) {
    vector h[];
    heights(f,h);
    scalar pos[];
    position (f, pos, {0,1});
    double max = statsf(pos).max;
    double min = statsf(pos).min;
	static FILE * fp0 = fopen("pressure/waveprobe0.dat", "w");
	double ycoords[2]  = {0,max}; // defines the minimum and maximum y-value 
	double yMax0 = wprobe(4*r_0,ycoords,100); 
	// update file
	fprintf(fp0, "%g %g\n",t,yMax0);
	fflush(fp0);
}


event droplets (t += 0.1)
{
    scalar m[];
    foreach()
        m[] = f[] > 1e-3;
    int n = tag (m);

    double v[n];
    coord b[n];
    for (int j = 0; j < n; j++)
        v[j] = b[j].x = b[j].y = b[j].z = 0.;
    foreach_leaf()
        if (m[] > 0) {
        int j = m[] - 1;
        v[j] += dv()*f[];
        coord p = {x,y,z};
        foreach_dimension()
        b[j].x += dv()*f[]*p.x;
        }

    #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif

    char name[80];
    sprintf (name, "interface/bubbles");
    static FILE * fp = fopen (name, "w");

    for (int j = 0; j < n; j++)
        fprintf (fp, "%d %g %d %g %g %g\n", i, t,
            j, v[j], b[j].x/v[j], b[j].y/v[j]);
    fflush (fp);
}
*/

