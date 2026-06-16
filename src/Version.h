#ifndef VERSION_H
#define VERSION_H
#include "Gordion.h"

/*
Version history:

1.0 Jan 28, 2026
    Initial release

1.1 Feb 20, 2026
    Added version numbers to save file headers.
    Added VEV32 preprocessor symbol & real typedef
	for selecting 32 bit vev's and poly coeffs.

1.2 Apr 11, 2026
    Fixed Version::arith and added Version::rsize.
    Removed installation of even-length fermion generators in Gen::geninit.
    Changed autoEgens to autoToddgens, and added Op::flipT to auto-create
	T-odd fermion generators from T-even generator definitions.
    Added Numerics::tikhanov to allow adjustment of Lagrange matrix regularization.

1.3 May 24, 2026
    Fixed bogus Gauss law use in 1D Efermion canonicalizations.
    Changed spectrum eigensystem to use smaller second order form.
    Removed Numerics::tikhanov shift as not useful.
    Added 'active' boolean to Gen's.
    Added generator add, suspend & activate commands.
    Deprecated previous "add generator" command.
    Support use for Eigen library instead of Armadillo for linear algebra,
	intended just for testing with 80-bit precision on Intel processors.
    Fixed Global::stageinit and related routines so that observable lists and
	polynomials scripts are not unnecessarily cleared upon a change of stage.

1.4 June 6, 2026
    Revamped initialization routines to disentangle initialization and generator
	definition from minimization stage. Changes include creating the OpList class
	and splitting Op::list into separate gauge and fermion OpList's held in
	global.stageinfo.ops.
    Moved nobsG and nobsF members from Global into associated ObsList.
    Ensure that Obs::catalog does not store into a frozen list.
    Added "stage" member to Coupling class, and made setting of coupling value
	automatically switch to associated minimization stage.
    Fixed initial setting of fermion maxord to 1 in obsinit, else GAxf misclassified.
    Fixed critical error in spec obs addition to Canon::cache.
*/

class Version					// Program version
    {
    static constexpr uchar	majornum = 1 ;
    static constexpr uchar	minornum = 4 ;

    static constexpr bool	r_is_f { sizeof (real) == sizeof (float) } ;
    static constexpr ushort	mkversion { (majornum << 4) | (minornum << 1) | r_is_f } ;

    ushort	version = mkversion ;

    ushort	rsize () const { return sizeof (real) ; } ;
    ushort	major () const { return (version >> 4) ; }
    ushort	minor () const { return (version >> 1) & 7 ; }
    bool	arith () const { return (version &  1) ; }

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
