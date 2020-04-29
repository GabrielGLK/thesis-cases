/**
# Momentum-conserving formulation for two-phase interfacial flows

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid
2. The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively. */

#include "all-mach-pc.h"
#include "vof-pc.h"

scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ and average viscosity $\mu$ (on faces) as well
as the cell-centered density. */

face vector alphav[], muv[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  mu = muv;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(f)  (clamp(f,0,1)*(mu1 - mu2) + mu2)
#endif

event properties (i++) {
  // fixme: metric
  foreach()
    rhov[] = rho(f[]);
  boundary ({rhov});
  foreach_face () {
    alphav.x[] = 2./(rhov[] + rhov[-1]);
    double ff = (f[] + f[-1])/2.;
    muv.x[] = fm.x[]*mu(ff);
  }
  boundary ((scalar *){muv});
}

/**
We overload the *vof()* event to transport consistently the volume
fraction and the momentum of each phase. */

static scalar * interfaces1 = NULL;

event vof (i++) {

  /**
  We split the total momentum $q$ into its two components $q1$ and
  $q2$ associated with $f$ and $1 - f$ respectively. */
  
  vector q1 = q, q2[];
  foreach()
    foreach_dimension() {
      double u = q.x[]/rho(f[]);
      q1.x[] = f[]*rho1*u;
      q2.x[] = (1. - f[])*rho2*u;
    }
  boundary ((scalar *){q1,q2});

  /**
  Momentum $q2$ is associated with $1 - f$, so we set the *inverse*
  attribute to *true*. We use (strict) minmod slope limiting for all
  components. */

  theta = 1.;
  foreach_dimension() {
    q2.x.inverse = true;
    q1.x.gradient = q2.x.gradient = minmod2;
  }

  /**
  We associate the transport of $q1$ and $q2$ with $f$ and transport
  all fields consistently using the VOF scheme. */

  f.tracers = (scalar *){q1,q2};
  vof_advection ({f}, i);

  /**
  We recover the total momentum. */
  
  foreach()
    foreach_dimension()
      q.x[] = q1.x[] + q2.x[];
  boundary ((scalar *){q});

  /**
  We set the list of interfaces to NULL so that the default *vof()*
  event does nothing (otherwise we would transport $f$ twice). */
  
  interfaces1 = interfaces, interfaces = NULL;
}

/**
We set the list of interfaces back to its default value. */

event tracer_advection (i++) {
  interfaces = interfaces1;
}

/**
## See also

* [Two-phase interfacial flows](two-phase.h)
*/
