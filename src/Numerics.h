#ifndef NUMERICS_H
#define NUMERICS_H
#include "Global.h"
#include "Poly.h"
#include "Ode.h"

using Numvec = arma::Col<real> ;
using Nummtx = arma::mat ;
using Cmplxv = arma::cx_vec ;
using Cmplxm = arma::cx_mat ;
using uvec   = arma::uvec ;

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
    Numvec	gradient ;		// Hamiltonian gradient
    Nummtx	curvature ;		// Hamiltonian curvature
    Numvec	delta ;			// Riemann normal coords
    Numvec	vev ;			// Numerical expectation values
    Numvec	dvev ;			// Expectation value derivatives
    Nummtx	lagrange ;		// Lagrange matrix
    Cmplxv	spectrum ;		// Oscillation spectrum
    Cmplxm	modes ;			// Oscillation eigenvectors
    short	lastrep ;		// Symmetry representation

    Numvec	vev_tmp ;		// Temporary vev vector
    const doub*	vev_buf ;		// Pointer to vev buffer data
    doub*	dvev_buf ;		// Pointer to dvev buffer data

    doub	dflttol = 1.e-8 ;	// Default tolerance
    doub	dfltlim = 1.e-7 ;	// Default svd limit
    doub	mintol  = dflttol ;	// Minimiization tolerance
    doub	odetol  = dflttol ;	// ODE integration tolerance
    doub	svdlim  = dfltlim ;	// Singular value threshold
    doub	tikhonov = 0 ;		// Lagrange Tikhonov shift
    uint	odemax	= Ode::dfltmax ;// Max ODE integration steps
    uint	minmax  = 500 ;		// Max Newton iterations
    short	minlim  = 0 ;		// Minimization generator order limit
    short	speclim = 0 ;		// Spectrum generator order limit
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
    Numvec&	eval_grad	(bool=false) ;		// Evaluate dH
    Numvec&	eval_delta	(bool=false) ;		// Predict minimum
    Cmplxv&	eval_spectra	(int,bool=false) ;	// Evaluate spectrum
    Nummtx&	eval_lagr	(int,bool=false) ;	// Evaluate Lagrange brkt
    Nummtx&	eval_curv	(int,bool,int=0) ;	// Evaluate ddH
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
    static void	numericsinit () ;				// Initialize
    static bool	curv_rpt     (int = -1) ;			// Curvature probs?
    static void status_rpt   (uint,uint) ;			// Report status
    static void	do_dvev      (doub, const Numvec&, Numvec&) ;	// Do vev derivs
    static void	do_dvev_bckt (const uint3&) ;			// Do dvev bucket
    static bool check_curv   (const Nummtx&) ;			// Curvature OK?
    static doub	err_norm     (const Numvec&, const Numvec&) ;	// ODE error norm
    static bool built_rep    (int) ;				// Is rep built?
    static string MMAform    (doub) ;				// "E" -> "*^"

    static void	data_write   (ofstream&, const string, doub, doub) ;
    static void	data_write   (ofstream&, const string, doub, const Numvec&) ;
    static void	data_write   (ofstream&, const string, doub, const Cmplxv&) ;
    } ;

extern Numerics numerics ;

#endif
