#ifndef RK_H
#define RK_H
#include <limits>
#include <cstring>
#include <array>
#include "Gordion.h"				// just for real typedef

class RKdef					// embedded Runge-Kutta method
    {
    public:
    int			nstage ;	// # stages
    int			order ;		// Method order
    const real*		c ;		// RK nodes
    const real* const	*a ;		// RK matrix
    const real*		b ;		// RK weights (higher order)
    const real*		bh ;		// RK weights (lower order)
    const real*		db ;		// RK error weights (bh - b)
    const char*		name ;		// Method name
    bool		fsal ;		// First deriv same as last?

    RKdef (int n, int o,		// Constructor
	   const real* c,  const real* const* a,
	   const real* b,  const real* bh,
	   const real* db, const char *nam)
	:
	    nstage(n), order(o),
	    c(c), a(a), b(b), bh(bh), db(db),
	    name(nam), 
	    fsal {b[n-1] == 0 && !std::memcmp(a[n-2],b,(n-1)*sizeof(real))}
	{}

    static const std::array<RKdef,7>	list ;	// Array of defined methods
    } ;

#endif
