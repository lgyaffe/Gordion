#ifndef RK_H
#define RK_H
#include <limits>
#include <cstring>
#include <array>

class RKdef					// embedded Runge-Kutta method
    {
    public:
    int				nstage ;	// # stages
    int				order ;		// Method order
    const double*		c ;		// RK nodes
    const double* const		*a ;		// RK matrix
    const double*		b ;		// RK weights (higher order)
    const double*		bh ;		// RK weights (lower order)
    const double*		db ;		// RK error weights (bh - b)
    const char*			name ;		// Method name
    bool			fsal ;		// First deriv same as last?

    RKdef (int n, int o,			// Constructor
	   const double* c,  const double* const* a,
	   const double* b,  const double* bh,
	   const double* db, const char *nam)
	:
	    nstage(n), order(o),
	    c(c), a(a), b(b), bh(bh), db(db),
	    name(nam), 
	    fsal {b[n-1] == 0 && !std::memcmp(a[n-2],b,(n-1)*sizeof(double))}
	{}

    static const std::array<RKdef,7>	list ;	// Array of defined methods
    } ;

#endif
