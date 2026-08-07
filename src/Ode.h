#ifndef ODE_H
#define ODE_H
#include "Gordion.h"
#include "Linalg.h"
#include "RK.h"

using OdeRhs  = void (*)(doub, const Rvec&, Rvec&) ;
using OdeNorm = doub (*)(const Rvec&, const Rvec&) ;

class Ode : public RKdef			// Runge-Kutta ODE integrator
    {
    OdeRhs		rhs ;			// derivative function
    OdeNorm		norm ;			// error norm function

    public:
    doub		tolerance ;			// Error tolerance
    doub		minstepsize  = 0.0 ;		// Minimum step size
    doub		maxstepsize  = 0.0 ;		// Maximum step size
    doub		initstepsize = 0.0 ;		// Initial step size
    uint		steps        = 0 ;		// # steps taken
    uint		rejects      = 0 ;		// # rejceted steps
    uint		maxstep ;			// Max # steps

    bool		integrate (doub&, doub, Rvec&) ;	// Integrator

    Ode (OdeRhs rhsfunc, OdeNorm normfunc, doub& tol, const RKdef& rk, uint max=dfltmax)
	: rhs(rhsfunc), norm(normfunc), tolerance(tol), RKdef(rk), maxstep(max) {} ;

    static struct Stats
	{
	static inline uint steps   { 0 } ;		// # total steps
	static inline uint rejects { 0 } ;		// # steps rejected
	static inline uint derivs  { 0 } ;		// # derivative evals
	static inline uint integs  { 0 } ;		// # integrations
	static inline uint fails   { 0 } ;		// # failed integrations
	} stats ;

    static constexpr int dfltmax = 1500 ;	// Default max steps
    } ;

#endif
