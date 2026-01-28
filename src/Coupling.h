#ifndef COUPLING_H
#define COUPLING_H
#include "Obs.h"

class Coupling : public char8		// Adjustable coupling constant
    {
    public:
    doub	value = 0 ;			// Coupling value

    static int indx (const string& s)		// Return coupling index
	{
	if (s.size() < sizeof (char8))
	    for (int i(0) ; i < list.size() ; ++i)
		if (s == list[i].data()) return i ;
	return -1 ;
	}

    static int indx (const char8& name)		// Return coupling index
	{
	auto ptr { std::find (list.begin(), list.end(), name) } ;
	return ptr != list.end() ? ptr - list.begin() : -1 ;
	}

    Coupling (char8 name) : char8(name) {}	// Constructor
    Coupling () {}

    static vector<Coupling> list ;		// Coupling list

    friend ostream& operator<< (ostream& stream, const Coupling& coup)
	{ return stream << coup.data() ; }
    } ;

using Couplings = vector<Coupling> ;

class AdjTerm                           // ObsPoly times adjustable coupling
    {
    public:
    short	coupindx ;			// Coupling index
    short	exponent = 1 ;			// Coupling exponent
    bool	imag  ;				// Imaginary coefficient?
    ObsPoly	poly  {ObsList::base} ;		// Observable polynomial
    ObsPoly	cpoly {ObsList::obs} ;		// Canonicalized form

    AdjTerm (int i, int e, ObsPoly& p, ObsPoly &c, bool img = false)
	: coupindx(i), exponent(e), poly(p), cpoly(c), imag(img) {}

    friend ostream& operator<< (ostream&, const vector<AdjTerm>&) ;
    } ;

#endif
