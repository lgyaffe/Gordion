#ifndef TERM_H
#define TERM_H
#include "Print.h"

template <typename C,typename T> class Term 	// T object w. C coefficient
    {
    public:
    C	coeff = 0.0 ;
    T	item ;

    Term () : item() {}
    Term (const T t, C d = 1.0) : item(t), coeff(d) {}

    friend ostream& operator<< (ostream& stream, const Term<C,T>& x)
	{
	Print::coeffprt (stream, x.coeff) ;
	return  stream << x.item ;
	}
    } ;

#endif
