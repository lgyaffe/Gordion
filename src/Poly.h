#ifndef POLY_H
#define POLY_H
#include "Gordion.h"
#include "Data.h"
#include "Term.h"
#include "Theory.h"
#include "Gripe.h"
#include <functional>
#include <cstdlib>

class Gen ;
class Obs ;
class ObsList ;
class ObsPoly ;

class Polyindx : public array<uint,PSIZ>		// Obs index tuples
    {
    public:
    Polyindx (uint i = 0)     : array<uint,PSIZ>{i} {}
    Polyindx (uint i, uint j) : array<uint,PSIZ>{i,j} { mysort() ; }

    void mysort()					// sort indices
	{ std::sort (begin(), end(), std::greater()) ; }

    bool validate() const				// check ordering
	{
	int i (PSIZ) ;
	while (--i > 0 && (*this)[i-1] >= (*this)[i]) ;
	return i == 0 ;
	}
    int order() const					// monomial order
	{
	for (int i (PSIZ-1) ; i >= 0 ; --i)
	    if ((*this)[i]) return i+1 ;
	return 0 ;
	}
    friend ostream& operator<< (ostream&, const Polyindx&) ;
    } ;

class PolyTerm : public Term<real,Polyindx>		// Polynomial term
    {
    public:
    using	Term<real,Polyindx>::Term ;

    operator	Polyindx()	  const { return item ; }	// Conversion
    const uint	operator[](int i) const { return item[i] ; }	// Subscript
    uint&	operator[](int i)	{ return item[i] ; }	// Subscript

    PolyTerm	operator* (const PolyTerm&) const ;		// Combine
    PolyTerm	operator* (doub z) const			// Scale term
		{ PolyTerm ans { *this } ; ans.coeff *= z ; return ans ; }

    int		order()    const { return item.order() ; }	// Term order
    bool	validate() const { return item.validate() ; }
    void	print (ostream&, const ObsList&) const ;
    } ;

struct Polyhash						// Polyindx hash function
    {
    std::size_t operator() (const Polyindx&) const noexcept ;
    } ;

class PolyMap : public hash<Polyindx,doub,Polyhash>	// PolyTerm hash table
    {
    std::reference_wrapper<ObsList> list ;		// Observable list

    public:
    PolyMap (ObsList& l) : list(l) {}			// Constructor

    ObsList& obslist() const { return list ; }		// Underlying ObsList
    
    void add (Polyindx t, doub d)			// Add term to PolyMap
	{
	auto [iter, isnew] { try_emplace (t, d) } ;
	if (!isnew) iter->second += d ;
	}
    void add (const PolyTerm&& t) { add (t.item, t.coeff) ; }	// Add term

    bool add_gen (const Gen&) ;				// Add Gen to PolyMap

    PolyMap& negate ()					// Negate entries
	{
	for (auto& entry : *this) entry.second *= -1 ;
	return *this ;
	}
    void purge ()					// Purge zero entries
	{
	std::erase_if (*this, [](const pair<Polyindx,doub>& p)
	    { return !p.second ; }) ;
	}
    bool allzero () { purge() ; return !size() ; }	// Vanishing Poly?

    friend ostream& operator<< (ostream&, const PolyMap&) ;
    } ;

class ObsPoly : public vector<PolyTerm>			// Cubic polynomial of Obs
    {
    public:
    std::reference_wrapper<ObsList> list ;		// Observable list

    void	push_map (PolyMap&) ;			// Add terms in map
    void	sort() ;				// Sort terms
    bool	allzero() const ;			// Vanishnig poly?
    ObsPoly&	scale(doub) ;				// Scale poly
    ObsPoly&	negate() { return scale(-1.0) ; }	// Negate poly
    ObsList&	obslist() const { return list ; }	// Underlying ObsList

    ObsPoly (const Gen&, ObsList&) ;			// Constructor
    ObsPoly (Obs&, ObsList&) ;				// Constructor
    ObsPoly (ObsList& l) : list(l) {}			// Constructor
    ObsPoly (uint indx, ObsList& l) : list(l),		// Constructor
		vector<PolyTerm>(1, PolyTerm(indx)) {}

    friend ostream& operator<< (ostream&, const ObsPoly&) ;
    } ;

struct Poly : public RecHdr				// Packed ObsPoly
    {
    public:
    const PolyTerm* begin() const
	{ return cast_to<const PolyTerm*>(this + 1) ; }

    const PolyTerm* end() const
	{ return cast_to<const PolyTerm*>(this + 1 + len * ptermsize) ; }

    friend ostream& operator<< (ostream&, const Poly&) ;

    static constexpr int ptermsize = sizeof (PolyTerm) / sizeof (RecHdr) ;
    } ;

static_assert (sizeof (PolyTerm) % sizeof (Poly) == 0) ;

class PolyRec : public DataRec				// Polynomial data record
    {
    public:
    mutable vector<ulong> offset ;			// Item offests

    using DataRec::DataRec ;
    PolyRec (DataRec rec) : DataRec (rec) {}		// Copy constructor

    void add (RecHdr, PolyMap&) ;				// Add PolyTerms

    void clear () { DataRec::clear() ; offset.clear() ; }	// Clear data

    const Poly& operator() (uint k, uint j=0, uint i=0) const	// Indexed Poly
	{						// N.B. slice/col major!
	if (offset.empty()) slice_n_dice () ;
	auto dataptr { begin().ptr } ;
	return dataptr[offset [(i * entry().ncol + j) * entry().nrow + k]] ;
	}

    void slice_n_dice () const				// Initialize offset array
	{
	auto n { entry().items() } ;

	offset.reserve (n) ;
	for (const auto& poly : *this)
	    {
	    offset.push_back (&poly - begin().ptr) ; --n ;
	    }
	if (n) gripe ("Corrupted data block!") ;
	}

    struct PolyIter					// Packed Poly iterator
	{
	const Poly*	ptr ;
	const Poly&	operator*() { return *ptr ; }
	PolyIter& 	operator++()
			{ ptr += 1 + ptr->len * Poly::ptermsize ; return *this ; }
	bool		operator!= (const PolyIter& a) { return ptr != a.ptr ; }

	PolyIter (const Poly* p) : ptr(p) {}
	} ;

    const PolyIter begin() const
	{ return PolyIter (static_cast<const Poly*>(data())) ; }

    const PolyIter end()   const
	{ return PolyIter (static_cast<const Poly*>(data()) + size()) ; }

    } ;

template <size_t N>
class PolyArr : public array<PolyRec,N>			// Array of PolyRec's
    { public: PolyArr (RecordID id) ; } ;

#endif
