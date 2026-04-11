#include "Ode.h"
#include <cmath>

bool Ode::integrate (double& s, double goal, Numvec& y)	//  Integrate dy/ds = rhs(y,s)
    {
    static Numvec		trial ;
    static Numvec		error ;
    static std::vector<Numvec>	deriv ;

    if (deriv.size() != nstage)   deriv.resize (nstage) ;
    if (trial.size() != y.size()) trial.set_size (y.size()) ;
    if (error.size() != y.size()) error.set_size (y.size()) ;
    for (auto& b : deriv)
	if (b.size() != y.size()) b.set_size (y.size()) ;

    if (!maxstepsize || maxstepsize > abs(goal - s)) maxstepsize = abs (goal - s) ;

    double	stepsize { initstepsize ? initstepsize : maxstepsize } ;
    double	errnorm   (0) ;
    int		stepcount (0) ;
    int		failures  (0) ;
    bool	known_deriv { false } ;

    while (s != goal)
	{
	if (++stepcount > maxstep) break ;

	if (errnorm) stepsize *= (failures > 1) ? 0.5 : (errnorm > tolerance) ?
				 std::max (0.9 * pow (tolerance/errnorm, 1./order),     0.5) :
				 std::min (0.9 * pow (tolerance/errnorm, 1./(order-1)), 2.0) ;
	if (minstepsize && stepsize < minstepsize) stepsize = minstepsize ;
	if (maxstepsize && stepsize > maxstepsize) stepsize = maxstepsize ;

	double next_s  { goal } ;
	double delta_s { abs (next_s - s) } ;
	if (stepsize < delta_s)
	    {
	    if (stepsize > delta_s / 2) stepsize = delta_s / 2 ;
	    next_s = (goal > s) ? s + stepsize : s - stepsize ;
	    }
	double h { next_s - s } ;

	if (!known_deriv) { rhs (s, y, deriv[0]) ; ++stats.derivs ; }

	for (int i (1) ; i < nstage ; ++i)
	    {
	    trial = y ;
	    for (int j(0) ; j < i ; ++j)
		if (a[i-1][j]) trial += h * a[i-1][j] * deriv[j] ;

	    rhs (s + c[i] * h, trial, deriv[i]) ;
	    ++stats.derivs ;
	    }
	if (!fsal)
	    {
	    trial = y ;
	    for (int i(0) ; i < nstage ; ++i)
		if (b[i]) trial += h * b[i] * deriv[i] ;
	    }
	error.zeros() ;
	for (int i(0) ; i < nstage ; ++i)
	    if (db[i]) error += h * db[i] * deriv[i] ;

	errnorm = norm (error, trial) ;
	if (errnorm < tolerance)
	    {
	    s = next_s ;
	    y = trial ;
	    failures = 0 ;
	    if (fsal)
		{
		deriv.front() = deriv.back() ;
		known_deriv = true ;
		}
	    else known_deriv = false ;
	    }
	else
	    {
	    known_deriv = true ;
	    ++failures ;
	    ++rejects  ;
	    ++stats.rejects ;
	    }
	++steps ;
	++stats.steps ;
	}
    ++stats.integs ;
    if (s == goal) return true ;
    ++stats.fails ;
    return false ;
    }
