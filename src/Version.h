#ifndef VERSION_H
#define VERSION_H
#include "Gordion.h"

/*
Version history:

1.0 Jan 28, 2026 initial release
1.1 Feb 20, 2026 added version numbers to save file headers,
		 added VEV32 preprocessor symbol & real typedef
		 for selecting 32 bit vev's and poly coeffs
*/

class Version					// Program version
    {
    static constexpr uchar	majornum = 1 ;
    static constexpr uchar	minornum = 1 ;

    static constexpr bool	r_is_f { sizeof (real) == sizeof (float) } ;
    static constexpr ushort	myversion { (majornum << 4) | (minornum << 1) | r_is_f } ;

    ushort	major () const { return (version >> 4) ; }
    ushort	minor () const { return (version >> 1) & 7 ; }
    ushort	arith () const { return (version &  1) ? sizeof (float) : sizeof (doub) ; }

    ushort	version = myversion ;

    public:
    bool incompat () const { return arith() != r_is_f ; }
    bool newer	  () const { return (version >> 1) > (myversion >> 1) ; }

    friend ostream& operator<< (ostream& stream, const Version v)
	{
	stream << v.major() << "." << v.minor() << "." << 8 * v.arith() ;
	return stream ;
	}
    } ;

#endif
