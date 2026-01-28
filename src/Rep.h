#ifndef REP_H
#define REP_H
#include "Symm.h"
#include "Gen.h"

class Proj : public vector<doub>		// Symmetry irrep projector
    {
    public:
    string		name ;			// representation name
    doub		denom = 1 ;		// normalization

    Proj () : vector<doub> (Symm::list.size(), 0) {}	// Constructor
    Proj (const string&, doub, SymmSum&&) ;		// Constructor

    Proj	operator()(const Proj&) const ;	// Compose projectors
    Proj	operator+ (const Proj&) const ;	// Add projectors
    Proj	operator- (const Proj&) const ;	// Subtract projectors
    Proj&	operator+=(const Proj&) ;	// Add to projector
    Proj&	operator-=(const Proj&) ;	// Subtract from projector
    bool	allzero() const ;		// Test for all zeros
    bool	C_even()  const ;		// Charge conjugation even?
    uint2	indices () const ;		// Projector indices
    strview	rowname () const ;		// Rep row name
    strview	repname () const ;		// Rep name

    static vector<Proj> list ;			// Known projectors

    friend ostream& operator<< (ostream&, const Proj&) ;
    } ;

class Rep : public vector<Proj>			// Representation projector(s)
    {
    public:
    using vector<Proj>::vector ;
    string	name ;				// representation name
    bool	C_even ;			// Charge conjugation even?

    Rep (string n, bool C) : name(n), C_even(C) {}

    static Index<Rep>	list ;			// Known representations
    static uint		known(const string&) ;	// Known Rep?
    static void		repinit() ;		// Initialize

    friend ostream& operator<< (ostream&, const Rep&) ;
    } ;

#endif
