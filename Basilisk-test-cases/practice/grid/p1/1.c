#include "grid/cartesian.h"
scalar f[];

/*
f[left] = dirichlet(0.);
f[right] = dirichlet(0.);
f[top] = dirichlet(0.);
f[bottom] = dirichlet(0.);
*/
int main()
{
    size(1);
    init_grid(10);
    foreach()
        f[] = sin(2.*pi*x) + cos(2.*pi*y);
    boundary({f});

    scalar dfx1[],dfx2[],dfx3[],dfx4[];
    foreach()
        {
            dfx1[] = (f[1] - f[])/Delta;//Forward difference scheme
            dfx2[] = (f[] - f[-1])/Delta;//Backward difference scheme
            dfx3[] = (f[1] - f[-1])/(2*Delta);//centered difference scheme
            dfx4[] = (f[1] - 2*f[] + f[-1])/(sq(Delta));
        }
    boundary({dfx1,dfx2,dfx3,dfx4});

    scalar dfy[];
    foreach()
        dfy[] = (f[0,1] - f[])/Delta;
    boundary({dfy});

    face vector uf[];
    foreach_face()
        uf.x[] = (f[1] - f[])/Delta;
    boundary((scalar *){uf});

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 1024);
    sprintf (name, "log");
    FILE *fp1 = fopen (name, "w");
    foreach()
        fprintf(fp1,"%g %g %g %g %g %g %g\n",x,y,f[],dfx1[],dfx2[],dfx3[],dfx4[]);

    sprintf (name, "log1");
    FILE *fp2 = fopen (name, "w");
    foreach_face()
        fprintf(fp2,"%g %g %g %g\n",x,y,uf.x[], uf.y[]); // note uf.y[] and uf.x[] are two inverse face vector

    

    name = (char *) malloc(sizeof(char) * 1024);
    sprintf (name, "cell");
    FILE *fp = fopen (name, "w");
    output_cells (fp);
    
    free_grid();
}