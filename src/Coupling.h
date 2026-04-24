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

struct Factor
    {
    int		indx ;
    doub	exp ;
    } ;

class Coeff : public vector<Factor>
    {
    public:
    using vector::vector ;
    doub operator() () const
	{
	auto&	list { Coupling::list } ;
	doub	ans  { 1.0 } ;
	for (auto& f : *this)
	    if (f.exp) ans *= std::pow(list[f.indx].value, f.exp) ;
	return ans ;
	}
    } ;

class AdjTerm                           // ObsPoly times adjustable coupling
    {
    public:
    Coeff	coeff ;				// Adjustable coefficient
    bool	imag  ;				// Imaginary coefficient?
    ObsPoly	poly  {ObsList::base} ;		// Observable polynomial
    ObsPoly	cpoly {ObsList::obs} ;		// Canonicalized form

    AdjTerm (Coeff& c, ObsPoly& p, ObsPoly &can, bool img = false)
	: coeff(c), poly(p), cpoly(can), imag(img) {}

    friend ostream& operator<< (ostream&, const vector<AdjTerm>&) ;
    } ;

#endif
