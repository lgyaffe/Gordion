#ifndef ODE_H
#define ODE_H
#define ARMA_WARN_LEVEL 1
#include <armadillo>
#include "RK.h"

using Numvec  = arma::vec ;
using OdeRhs  = void   (*)(double, const Numvec&, Numvec&) ;
using OdeNorm = double (*)(const Numvec&, const Numvec&) ;

class Ode : public RKdef			// Runge-Kutta ODE integrator
    {
    OdeRhs		rhs ;			// derivative function
    OdeNorm		norm ;			// error norm function

    public:
    double		tolerance ;			// Error tolerance
    double		minstepsize  = 0.0 ;		// Minimum step size
    double		maxstepsize  = 0.0 ;		// Maximum step size
    double		initstepsize = 0.0 ;		// Initial step size
    uint		steps        = 0 ;		// # steps taken
    uint		rejects      = 0 ;		// # rejceted steps
    uint		maxstep ;			// Max # steps

    bool		integrate (double&, double, Numvec&) ;	// Integrator

    Ode (OdeRhs rhsfunc, OdeNorm normfunc, double& tol, const RKdef& rk, uint max=dfltmax)
	: rhs(rhsfunc), norm(normfunc), tolerance(tol), RKdef(rk), maxstep(max) {} ;

    struct Stats
	{
	static inline uint32_t steps   { 0 } ;		// # total steps
	static inline uint32_t rejects { 0 } ;		// # steps rejected
	static inline uint32_t derivs  { 0 } ;		// # derivative evals
	static inline uint32_t integs  { 0 } ;		// # integrations
	static inline uint32_t fails   { 0 } ;		// # failed integrations
	} ;

    static constexpr int	dfltmax = 1500 ;	// Default max steps
    static inline Stats		stats ;			// ODE integration stats
    } ;

#endif
