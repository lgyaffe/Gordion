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
    array<dcmplx,3>	cmplxeig ;	// Complex curvature eigs
    doub		maxloopv ;	// Max loop expectation
    numb		maxloopi ;	// Max loop index

    void reset () { warn = false ; negeig.fill(0) ; cmplxeig.fill(0) ; }
    } ;

class Numerics
    {
    public:
    Dvec	gradient ;		// Hamiltonian gradient
    Dmtx	curvature ;		// Hamiltonian curvature (T-even)
    Dmtx	metric ;		// Hamiltonian curvature (T-odd)
    Dmtx	lagrange ;		// Lagrange bracket matrix
    Dvec	delta ;			// Riemann normal coords
    Rvec	vev ;			// Numerical expectation values
    Rvec	dvev ;			// Expectation value derivatives
    Cdvec	spectrum ;		// Oscillation spectrum
    Cdmtx	modes ;			// Oscillation eigenvectors
    Uvec	Tevens ;		// T-even active generators
    Uvec	Todds ;			// T-odd active generators
    ushort	lastrep ;		// Symmetry representation

    uint	maxode	 { Ode::dfltmax } ; // Max ODE steps
    uint	maxnewt  { 100 } ;	// Max Newton iterations
    bool	oknegeig { false } ;	// Negative curvature OK?
    bool	symcurv  { true } ;	// Symmetrize curvature?
    RKdef	rk ;			// RK method

    Rvec	vev_tmp ;		// Temporary vev vector
    const real*	vev_buf ;		// Pointer to vev buffer data
    real*	dvev_buf ;		// Pointer to dvev buffer data

    struct
	{
	uint	tries = 0 ;		// # minimizations
	uint	fails = 0 ;		// # minimization failures
	}	stats ;			// Minimization stats

    void	 do_minimize	() ;			// Do H minimization
    int		 do_step	(doub=0) ;		// Do minimization step
    void	 do_flow	(uint,doub,doub,doub) ;	// Do coupling flow
    doub	 eval_ham	(bool=false) ;		// Evaluate H
    void	 eval_geos	(int=0) ;		// Evaluate Obs derivs
    doub	 eval_delta	(bool=false) ;		// Predict minimum
    const Dvec&	 eval_grad	(bool=false) ;		// Evaluate dH
    const Dmtx&	 eval_curv	(uint,int=0) ;		// Evaluate T-even ddH
    const Dmtx&	 eval_metr	(uint,int=0) ;		// Evaluate T-odd ddH
    const Dmtx&	 eval_lagr	(uint,int=0) ;		// Evaluate Lagrange brkt
    const Uvec&	 eval_inuse	(uint,bool=false) ;	// Active generator list
    const Cdvec& eval_spectra	(uint,bool=false) ;	// Evaluate spectrum
    const Cdvec& eval_spectra	(string,bool=false) ;	// Evaluate spectrum

    bool	check_loops	() ;			// Loop vevs < 1?
    void	initialize	(uint=global.stage) ;	// Initialize & reset
    void	status_rpt	(uint,uint) ;		// Report status
    bool	open_MMA	() ;			// Open MMA output file
    void	write_data	() ;			// Write to MMA file

    static doub termvalue (const PolyElem*& ptr, const real* v)	// Evaluate Poly
	{
	doub z	  { (ptr++)->coeff } ;
	numb indx { (ptr++)->index } ;
	for (; indx < 0 ; indx = (ptr++)->index) z *= v[-indx] ;
	return z * v[indx] ;
	}

    doub termvalue (const PolyElem*& ptr)			// Evaluate Poly
	 { return termvalue (ptr, vev.memptr()) ; }

#ifdef NUM64
    static constexpr doub DFLTTOL { 1.e-10 } ;	// Default ODE tolerance
#else
    static constexpr doub DFLTTOL { 1.e-6 } ;	// Default ODE tolerance
#endif
    static inline doub	 svdcut  { 0 } ;	// Singular value cutoff
    static inline doub	 dflttol { DFLTTOL } ;	// Default min & ODE tolerance
    static inline doub	 mintol  { dflttol } ;	// Minimiization tolerance
    static inline doub	 odetol  { dflttol } ;	// ODE integration tolerance
    static inline Status status ;		// Status info

    static bool		check_curv (const Dmtx&) ;	// Curvature OK?
    static bool		curv_rpt   (int = -1) ;		// Curvature probs?
    static bool 	built_rep  (uint) ;		// Is rep built?
    static string	MMAform    (doub) ;		// "E" -> "*^"

    static void	do_dvev      (doub, const Rvec&, Rvec&) ; // Do vev derivs
    static void	do_dvev_bckt (const numb3&) ;		  // Do dvev bucket
    static doub	err_norm     (const Rvec&, const Rvec&) ; // ODE error norm

    static void	data_write (ofstream&, const string, doub) ;
    static void	data_write (ofstream&, const string, const Rvec&) ;
    static void	data_write (ofstream&, const string, const Cdvec&) ;
    } ;

inline Numerics numerics ;				// Numerical data

#endif
