#ifndef VERSION_H
#define VERSION_H
#include "Gordion.h"

class Version					// Program version
    {
    static constexpr uchar	majornum = 1 ;
    static constexpr uchar	minornum = 5 ;

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
