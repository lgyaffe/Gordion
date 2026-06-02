#ifndef VERSION_H
#define VERSION_H
#include "Gordion.h"

/*
Version history:

1.0 Jan 28, 2026
    initial release

1.1 Feb 20, 2026
    added version numbers to save file headers.
    added VEV32 preprocessor symbol & real typedef
	for selecting 32 bit vev's and poly coeffs.

1.2 Apr 11, 2026
    fixed Version::arith and added Version::rsize.
    removed installation of even-length fermion generators in Gen::geninit.
    changed autoEgens to autoToddgens, and added Op::flipT to auto-create
	T-odd fermion generators from T-even generator definitions.
    added Numerics::tikhanov to allow adjustment of Lagrange matrix regularization.

1.3 May 24, 2026
    fixed bogus Gauss law use in 1D Efermion canonicalizations.
    changed spectrum eigensystem to use smaller second order form.
    removed Numerics::tikhanov shift as not useful.
    added 'active' boolean to Gen's.
    added generator add, suspend & activate commands.
    deprecated previous "add generator" command.
    support use for Eigen library instead of Armadillo for linear algebra.
    fixed Global::stageinit and related routines so that existing observable lists
    and polynomials scripts are not cleared upon a change of stage from gauge to fermi
    or back.  Polynomial scripts are only cleared if additional generators are defined.
    Fermion observables are only cleared if additional gauge observables are to be built.
*/

class Version					// Program version
    {
    static constexpr uchar	majornum = 1 ;
    static constexpr uchar	minornum = 3 ;

    static constexpr bool	r_is_f { sizeof (real) == sizeof (float) } ;
    static constexpr ushort	mkversion { (majornum << 4) | (minornum << 1) | r_is_f } ;

    ushort	rsize () const { return r_is_f ? sizeof (float) : sizeof (doub) ; } ;
    ushort	major () const { return (version >> 4) ; }
    ushort	minor () const { return (version >> 1) & 7 ; }
    bool	arith () const { return (version &  1) ; }

    ushort	version = mkversion ;

    public:
    bool incompat () const { return arith() != r_is_f ; }
    bool newer	  () const { return (version >> 1) > (mkversion >> 1) ; }

    friend ostream& operator<< (ostream& stream, const Version v)
	{
	stream << v.major() << "." << v.minor() << "." << 8 * v.rsize() ;
	return stream ;
	}
    } ;

#endif
