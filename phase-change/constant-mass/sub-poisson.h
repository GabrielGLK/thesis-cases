#include "curvature.h"

void mg_cycle_new (scalar * a, scalar * res, scalar * da,
	       void (* relax) (scalar * da, scalar * res, 
			       int depth, void * data),
	       void * data,
	       int nrelax, int minlevel, int maxlevel)
{

  /**
  We first define the residual on all levels. */

  restriction (res);

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel)
      foreach_level_or_leaf (l)
	for (scalar s in da)
	  s[] = 0.;

    /**
    On all other grids, we take as initial guess the approximate solution
    on the coarser grid bilinearly interpolated onto the current grid. */

    else
      foreach_level (l)
	for (scalar s in da)
	  s[] = bilinear (point, s);

    /**
    We then apply homogeneous boundary conditions and do several
    iterations of the relaxation function to refine the initial guess. */

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }

  /**
  And finally we apply the resulting correction to *a*. */

  foreach() { 
    scalar s, ds;
    for (s, ds in a, da)
    {if(interfacial(point,f)||f[] == 0)
      s[] += ds[];
    else
        s[] += 0;
    }
  }
  boundary (a);
}

/**
## Multigrid solver

The multigrid solver itself uses successive calls to the multigrid
cycle to refine an initial guess until a specified tolerance is
reached. 

The maximum number of iterations is controlled by *NITERMAX* and the
tolerance by *TOLERANCE* with the default values below. */


/**
Information about the convergence of the solver is returned in a structure. */

/**
The user needs to provide a function which computes the residual field
(and returns its maximum) as well as the relaxation function. The
user-defined pointer *data* can be used to pass arguments to these
functions. The optional number of relaxations is *nrelax* (default is
one) and *res* is an optional list of fields used to store the
residuals. The minimum level of the hierarchy can be set (default is
zero i.e. the root cell). */

struct MGSolve_New {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
		       void * data);
  void (* relax) (scalar * da, scalar * res, int depth, 
		  void * data);
  void * data;
  
  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats mg_solve_new (struct MGSolve_New p)
{

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    for (scalar s in p.a) {
      scalar r = new scalar;
      res = list_append (res, r);
    }

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
    {if(interfacial(point,f)||f[] == 0)
      sum += s[];
      else
        sum += 0;
    }
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */

  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data,
	      s.nrelax,
	      p.minlevel,
	      grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);

    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */

#if 1
    if (s.resa > TOLERANCE) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
	s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
	s.nrelax--;
    }
#else
    if (s.resa == resb) /* convergence has stopped!! */
      break;
    if (s.resa > resb/1.1 && p.minlevel < grid->maxdepth)
      p.minlevel++;
#endif

    resb = s.resa;
  }
  s.minlevel = p.minlevel;
  
  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr, 
	     "WARNING: convergence for %s not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d\n", v.name,
	     s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }
    
  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}

struct Poisson_New {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
#if EMBED
  double (* embed_flux) (Point, scalar, vector, double *);
#endif
};

/**
We can now write the relaxation function. We first recover the extra
parameters from the data pointer. */

static void relax_new (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;

  /**
  We use either Jacobi (under)relaxation or we directly reuse values
  as soon as they are updated. For Jacobi, we need to allocate space
  for the new field *c*. Jacobi is useful mostly as it gives results
  which are independent of the order in which the cells are
  traversed. This is not the case for the simple traversal, which
  means for example that results will depend on whether a tree or
  a multigrid is used (because cells will be traversed in a different
  order). The same comment applies to OpenMP or MPI parallelism. In
  practice however Jacobi convergence tends to be slower than simple
  reuse. */
  
#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif
  
  /**
  We use the face values of $\alpha$ to weight the gradients of the
  5-points Laplacian operator. We get the relaxation function. */

  foreach_level_or_leaf (l) {
    double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
    foreach_dimension() {
      n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
      d += alpha.x[1] + alpha.x[];
    }
#if EMBED
    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (!d)
      c[] = b[] = 0.;
    else
#endif // EMBED
      c[] = n/d;
  }

  /**
  For weighted Jacobi we under-relax with a weight of 2/3. */
  
#if JACOBI
  foreach_level_or_leaf (l)
    a[] = (a[] + 2.*c[])/3.;
#endif
  
#if TRASH
  scalar a1[];
  foreach_level_or_leaf (l)
    a1[] = a[];
  trash ({a});
  foreach_level_or_leaf (l)
    a[] = a1[];
#endif
}

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

static double residual_new (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*face_gradient_x (a, 0);
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
    {if(interfacial(point,f)||f[] == 0)
      res[] += (alpha.x[0]*face_gradient_x (a, 0) -
		alpha.x[1]*face_gradient_x (a, 1))/Delta; 
        else
            res[] += 0;}
#endif // !TREE    
#if EMBED
    if (p->embed_flux) {
      double c, e = p->embed_flux (point, a, alpha, &c);
      res[] += c - e*a[];
    }
#endif // EMBED    
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  boundary (resl);
  return maxres;
}

mgstats poisson_new (struct Poisson_New p)
{

  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction ({alpha,lambda});

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
#if EMBED
  if (!p.embed_flux && a.boundary[embed] != symmetry)
    p.embed_flux = embed_flux;
#endif // EMBED
  mgstats s = mg_solve_new ({a}, {b}, residual_new, relax_new,
			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
