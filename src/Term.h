#ifndef TERM_H
#define TERM_H
#include "Gordion.h"

extern ostream&	coeffprt (ostream&, doub) ;	// Print coefficient nicely

inline bool is_int (doub c)			// Is essentially integer?
    {
    return abs(c - std::round(c)) < 1.e-14 ;
    }

inline bool is_one (doub c)			// Is essentially 1?
    {
    return abs(c - 1) < 1.e-15 ;
    }

template <typename C,typename T> class Term 	// T object w. C coefficient
    {
    public:
    C	coeff = 0.0 ;
    T	item ;

    Term () : item() {}
    Term (const T t, C d = 1.0) : item(t), coeff(d) {}

    friend ostream& operator<< (ostream& stream, const Term<C,T>& x)
	{
	coeffprt (stream, x.coeff) ;
	return  stream << x.item ;
	}
    } ;

#endif
