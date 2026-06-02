#ifndef NUMERICS_H
#define NUMERICS_H
#include "Global.h"
#include "Poly.h"
#include "Ode.h"
#include "Linalg.h"

class Status
    {
    public:
    bool		warn ; 		// Warnings to report?
    array<doub,3>	negeig ;	// Negative curvature eigs
    array<cmplx,3>	cmplxeig ;	// Complex curvature eigs
    doub		maxloopv ;	// Max loop expectation
    uint		maxloopi ;	// Max loop index

    void reset () { warn = false ; negeig.fill(0) ; cmplxeig.fill(0) ; }
    } ;

class Numerics
    {
    public:
    doub	H ;			// Hamiltonian/free energy value
    Dvec	gradient ;		// Hamiltonian gradient
    Dmtx	curvature ;		// Hamiltonian curvature (T-even)
    Dmtx	metric ;		// Hamiltonian curvature (T-odd)
    Dmtx	lagrange ;		// Lagrange bracket matrix
    Dvec	delta ;			// Riemann normal coords
    Rvec	vev ;			// Numerical expectation values
    Rvec	dvev ;			// Expectation value derivatives
    Cvec	spectrum ;		// Oscillation spectrum
    Cmtx	modes ;			// Oscillation eigenvectors
    Uvec	Tevens ;		// T-even active generators
    Uvec	Todds ;			// T-odd active generators
    short	lastrep ;		// Symmetry representation

    Rvec	vev_tmp ;		// Temporary vev vector
    const real*	vev_buf ;		// Pointer to vev buffer data
    real*	dvev_buf ;		// Pointer to dvev buffer data

    doub	svdcut  = 0 ;		// Singular value cutoff
    doub	dflttol = 1.e-10 ;	// Default tolerance
    doub	mintol  = dflttol ;	// Minimiization tolerance
    doub	odetol  = dflttol ;	// ODE integration tolerance
    uint	odemax	= Ode::dfltmax ;// Max ODE integration steps
    uint	minmax  = 500 ;		// Max Newton iterations
    RKdef	rk {RKdef::list[0]} ;	// RK method

    struct
	{
	uint	tries = 0 ;		// # minimizations
	uint	fails = 0 ;		// # minimization failures
	}	stats ;			// Minimization stats

    void	do_minimize	() ;			// Do H minimization
    int		do_step		(doub=0) ;		// Do minimization step
    void	do_flow		(int,doub,doub,doub) ;	// Do coupling flow
    doub	eval_H		(bool=false) ;		// Evaluate H
    void	eval_geos	(int=0) ;		// Evaluate Obs derivs
    doub	eval_delta	(bool=false) ;		// Predict minimum
    const Dvec&	eval_grad	(bool=false) ;		// Evaluate dH
    const Dmtx&	eval_curv	(uint,int=0) ;		// Evaluate T-even ddH
    const Dmtx&	eval_metr	(uint,int=0) ;		// Evaluate T-odd ddH
    const Dmtx&	eval_lagr	(uint,int=0) ;		// Evaluate Lagrange brkt
    const Uvec&	eval_inuse	(uint,bool=false) ;	// Active generator list
    const Cvec&	eval_spectra	(uint,bool=false) ;	// Evaluate spectrum
    const Cvec&	eval_spectra	(string,bool=false) ;	// Evaluate spectrum
    void	write_data	(doub) ;		// Write to datafile

    inline static Status status ;			// Status info

    doub termvalue (const PolyTerm& t) const		// Evaluate PolyTerm
	{
	doub z { t.coeff } ;
	for (int k(0) ; k < PSIZ ; ++k)
	    if (t[k]) z *= vev[t[k]] ;
	    else break ;
	return z ;
	}

    static doub	termvalue (const PolyTerm& t, const doub* o)	// Evaluate PolyTerm
	{
	doub z { t.coeff } ;
	for (int k(0) ; k < PSIZ ; ++k)
	    if (t[k]) z *= o[t[k]] ;
	    else break ;
	return z ;
	}

    static bool check_loops  () ;				// Loop vevs < 1?
    static void	numericsinit (int = global.stage) ;		// Initialize
    static bool	curv_rpt     (int = -1) ;			// Curvature probs?
    static void status_rpt   (uint,uint) ;			// Report status
    static void	do_dvev      (doub, const Rvec&, Rvec&) ;	// Do vev derivs
    static void	do_dvev_bckt (const uint3&) ;			// Do dvev bucket
    static bool check_curv   (const Dmtx&) ;			// Curvature OK?
    static doub	err_norm     (const Rvec&, const Rvec&) ;	// ODE error norm
    static bool built_rep    (int) ;				// Is rep built?
    static string MMAform    (doub) ;				// "E" -> "*^"

    static void	data_write   (ofstream&, const string, doub, doub) ;
    static void	data_write   (ofstream&, const string, doub, const Rvec&) ;
    static void	data_write   (ofstream&, const string, doub, const Cvec&) ;
    } ;

extern Numerics numerics ;

#endif
