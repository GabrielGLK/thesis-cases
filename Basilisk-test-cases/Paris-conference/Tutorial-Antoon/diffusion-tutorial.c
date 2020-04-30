
#include "new-diffusion-2.h"
#include "run.h"

// we consider the evolution of a scalar filed s. It is decleared like so:
scalar s[];

/* we start the time loop by calling run().
we also set the value of DT, which will serve as a time-stepping parameter.
It became available by #include "run.h"
*/

#define LEVEL 8

int main() {
    init_grid(1 << (LEVEL - 2));
    DT = 0.005;
    run();
}

event init (t = 0){
    foreach()
        s[] = exp(-(sq(x-0.5) + sq(y-0.5))*10.);
    boundary({s});
}

event movie ( t += 0.1){
    output_ppm (s, file = "s.mp4", n = 256, min = -1, max = 1);
}


/*
we tell the time-loop to set the actual timestep(dt), based on our maximum value(DT). The
documentation for the diffusion solver reveals the proper sequence of the arguments to the 
diffusion() function.
*/
/* simple one
event diff (i++){
    dt = dtnext (DT);
    const face vector D[] = {0.01,0.01};
    diffusion( s, dt, D);
}*/

// avanced one
event diff( i++){
    dt = dtnext(DT);
    face vector D[];
    foreach_face(x)
        D.x[] = 0.01;
    foreach_face(y)
        D.y[] = 0.01;
    boundary((scalar *){D});// call boundary() function in order to ensure its proper definition near
    // resolutions boundaries and on the various lower levels.
    diffusion(s, dt, D);
}

#if 1 // The results are slightly different for uniform and quadtree grid
event adapt(i++)
    adapt_wavelet({s},(double[]){0.01}, 10);
#endif

event output(t += 1){
    char *name = NULL;
    name = (char *) malloc(sizeof(char) *1024);
    sprintf(name, "data-%g",t);
    static FILE *fp = fopen(name, "w");
    fprintf(fp, "%g %d %.8g\n", t, i, statsf(s).sum);
    //fclose(fp);  // this will cause "double free or corruption" error
}

event stop(t = 10);