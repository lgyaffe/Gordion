#ifndef NUMERICS_H
#define NUMERICS_H
#include "Global.h"
#include "Poly.h"
#include "Ode.h"

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
    ushort	lastrep ;		// Symmetry representation

    Rvec	vev_tmp ;		// Temporary vev vector
    const real*	vev_buf ;		// Pointer to vev buffer data
    real*	dvev_buf ;		// Pointer to dvev buffer data

    doub	svdcut  = 0 ;		// Singular value cutoff
    doub	dflttol = DFLTTOL ;	// Default min & ODE tolerance
    doub	mintol  = dflttol ;	// Minimiization tolerance
    doub	odetol  = dflttol ;	// ODE integration tolerance
    uint	maxode	= Ode::dfltmax ;// Max ODE integration steps
    uint	maxnewt = 500 ;		// Max Newton iterations
    RKdef	rk ;			// RK method

    struct
	{
	uint	tries = 0 ;		// # minimizations
	uint	fails = 0 ;		// # minimization failures
	}	stats ;			// Minimization stats

    void	do_minimize	() ;			// Do H minimization
    int		do_step		(doub=0) ;		// Do minimization step
    void	do_flow		(uint,doub,doub,doub) ;	// Do coupling flow
    doub	eval_ham	(bool=false) ;		// Evaluate H
    void	eval_geos	(int=0) ;		// Evaluate Obs derivs
    doub	eval_delta	(bool=false) ;		// Predict minimum
    const Dvec&	eval_grad	(bool=false) ;		// Evaluate dH
    const Dmtx&	eval_curv	(uint,int=0) ;		// Evaluate T-even ddH
    const Dmtx&	eval_metr	(uint,int=0) ;		// Evaluate T-odd ddH
    const Dmtx&	eval_lagr	(uint,int=0) ;		// Evaluate Lagrange brkt
    const Uvec&	eval_inuse	(uint,bool=false) ;	// Active generator list
    const Cvec&	eval_spectra	(uint,bool=false) ;	// Evaluate spectrum
    const Cvec&	eval_spectra	(string,bool=false) ;	// Evaluate spectrum
    void	initialize	(int = global.stage) ;	// Initialize
    void	status_rpt	(uint,uint) ;		// Report status
    bool	check_loops	() ;			// Loop vevs < 1?
    bool	open_MMA	() ;			// Open MMA output file
    void	write_data	() ;			// Write to MMA file

    doub termvalue (const PolyTerm& t)				// Evaluate PolyTerm
	 { return termvalue (t, vev.memptr()) ; }

    static doub	termvalue (const PolyTerm& t, const real* v)	// Evaluate PolyTerm
	{
	doub z { t.coeff } ;
	for (int k(0) ; k < PSIZ ; ++k)
	    if (t[k]) z *= v[t[k]] ;
	    else break ;
	return z ;
	}

    doub termvalue (const Poly*& ptr)				// Evaluate Poly
	 { return termvalue (ptr, vev.memptr()) ; }

    static doub termvalue (const Poly*& ptr, const real* v)	// Evaluate Poly
	{
	doub z    { (ptr++)->coeff } ;
	numb indx { (ptr++)->index } ;
	for (; indx < 0 ; indx = (ptr++)->index) z *= v[-indx] ;
	return z * v[indx] ;
	}

    inline static Status status ;				// Status info

    static bool check_curv   (const Dmtx&) ;			// Curvature OK?
    static bool	curv_rpt     (int = -1) ;			// Curvature probs?
    static bool built_rep    (uint) ;				// Is rep built?
    static string MMAform    (doub) ;				// "E" -> "*^"

    static void	do_dvev      (doub, const Rvec&, Rvec&) ;	// Do vev derivs
    static void	do_dvev_bckt (const numb3&) ;			// Do dvev bucket
    static doub	err_norm     (const Rvec&, const Rvec&) ;	// ODE error norm

    static void	data_write   (ofstream&, const string, doub) ;
    static void	data_write   (ofstream&, const string, const Rvec&) ;
    static void	data_write   (ofstream&, const string, const Cvec&) ;
    } ;

inline Numerics numerics ;				// Numerical data

#endif
