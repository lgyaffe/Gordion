#ifndef ODE_H
#define ODE_H
#include "Gordion.h"
#include "Linalg.h"
#include "RK.h"

using OdeRhs  = void (*)(doub, const Dvec&, Dvec&) ;
using OdeNorm = doub (*)(const Dvec&, const Dvec&) ;

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

    bool		integrate (doub&, doub, Dvec&) ;	// Integrator

    Ode (OdeRhs rhsfunc, OdeNorm normfunc, doub& tol, const RKdef& rk, uint max=dfltmax)
	: rhs(rhsfunc), norm(normfunc), tolerance(tol), RKdef(rk), maxstep(max) {} ;

    static struct Stats
	{
	static inline uint32_t steps   { 0 } ;		// # total steps
	static inline uint32_t rejects { 0 } ;		// # steps rejected
	static inline uint32_t derivs  { 0 } ;		// # derivative evals
	static inline uint32_t integs  { 0 } ;		// # integrations
	static inline uint32_t fails   { 0 } ;		// # failed integrations
	} stats ;

    static constexpr int dfltmax = 1500 ;	// Default max steps
    } ;

#endif
