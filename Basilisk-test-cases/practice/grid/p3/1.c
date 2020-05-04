#include "grid/cartesian1D.h"
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
        f[] = sin(2.*pi*x);
    boundary({f});

    scalar dfx1[],dfx2[],dfx3[],dfx4[],dfx5[],dfx6[],dfx7[];
    foreach()
        {
            dfx1[] = (f[1] - f[])/Delta;//Forward difference scheme
            dfx2[] = (f[] - f[-1])/Delta;//Backward difference scheme
            dfx3[] = (f[1] - f[-1])/(2*Delta);//centered difference scheme
            dfx4[] = (-3*f[] + 4*f[1] - f[2])/(2*Delta);
            dfx5[] = (4*f[1] + 6*f[] - 12*f[-1] + 2*f[-2])/(12*Delta);
            dfx6[] = (-2*f[2] + 12*f[1] - 6*f[] - 4*f[-1])/(12*Delta);
            dfx7[] = (f[-2] - 8*f[-1] + 8*f[1] - f[2])/(12*Delta);
            //dfx5[] = (-f[2] + 16*f[1] - 30*f[] + 16*f[-1] -f[-2])/(12*sq(Delta));//fourth order 2nd direvative differnce scheme
        }
    boundary({dfx1,dfx2,dfx3,dfx4,dfx5,dfx6,dfx7});

    char *name = NULL;
    name = (char *) malloc(sizeof(char) * 1024);
    sprintf (name, "log");
    FILE *fp1 = fopen (name, "w");
    foreach()
        fprintf(fp1,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",x,f[],dfx1[],dfx2[],dfx3[],dfx4[],dfx5[],dfx6[],dfx7[]);

    free_grid();
}