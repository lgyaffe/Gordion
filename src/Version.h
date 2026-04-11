#ifndef VERSION_H
#define VERSION_H
#include "Gordion.h"

/*
Version history:

1.0 Jan 28, 2026
    initial release

1.1 Feb 20, 2026
    added version numbers to save file headers,
    added VEV32 preprocessor symbol & real typedef
    for selecting 32 bit vev's and poly coeffs.

1.2 Apr 11, 2026
    fixed Version::arith and added Version::rsize,
    removed installation of even-length fermion generators in Gen::geninit,
    changed autoEgens to autoToddgens, and added Op::flipT to auto-create
    T-odd fermion generators from T-even generator definitions,
    added Numerics::tikhanov to allow adjustment of Lagrange matrix regularization.
*/

class Version					// Program version
    {
    static constexpr uchar	majornum = 1 ;
    static constexpr uchar	minornum = 2 ;

    static constexpr bool	r_is_f { sizeof (real) == sizeof (float) } ;
    static constexpr ushort	mkversion { (majornum << 4) | (minornum << 1) | r_is_f } ;

    ushort	rsize () const { return r_is_f ? sizeof (float) : sizeof (double) ; } ;
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
