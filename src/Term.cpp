#include "Term.h"
#include "Gordion.h"

ostream& coeffprt (ostream& stream, doub c)	// Nicely print object coefficient
    {
    if (c)
	{
	const string sgn { c > 0 ? "+" : "-" } ;
	c = abs(c) ;
	doub x { c * sqrt(2.0) } ;
	if      (is_one (c))	stream << sgn ;
	else if (is_int (c))	stream << sgn <<    c << " " ;
	else if (is_int (1/c))	stream << sgn << "1/" << (int)(1/c) << " " ;
	else if (is_int (2*c))	stream << sgn << 2*c  << "/2" << " " ;
	else if (is_int (3*c))	stream << sgn << 3*c  << "/3" << " " ;
	else if (is_int (4*c))	stream << sgn << 4*c  << "/4" << " " ;
	else if (is_int (8*c))	stream << sgn << 8*c  << "/8" << " " ;
	else if (is_int (16*c))	stream << sgn << 16*c << "/16" << " " ;
	else if (is_int (x))	stream << sgn <<   x  << "/\u221A2" << " "  ;
	else if (is_int (2*x))	stream << sgn << 2*x  << "/(2\u221A2)" << " " ;
	else if (is_int (3*x))	stream << sgn << 3*x  << "/(3\u221A2)" << " " ;
	else if (is_int (4*x))	stream << sgn << 4*x  << "/(4\u221A2)" << " " ;
	else if (is_int (8*x))	stream << sgn << 8*x  << "/(8\u221A2)" << " " ;
	else if (is_int (16*x))	stream << sgn << 16*x << "/(16\u221A2)" << " " ;
	else			stream << sgn << format("{:.2e} ", c) ;
	}
    else stream << "+0 " ;
    return stream ;
    } ;
