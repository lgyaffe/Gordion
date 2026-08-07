#ifndef POLY_H
#define POLY_H
#include "Data.h"
#include "Term.h"
#include "Theory.h"
#include <functional>
#include <utility>

class Gen ;
class ObsList ;
class PolyElem ;
class SysIndex ;

class PolyIndx : public array<numb,PSIZ>		// Obs index tuples
    {
    public:
    PolyIndx (numb i = 0)     : array<numb,PSIZ>{i} {}
    PolyIndx (numb i, numb j) : array<numb,PSIZ>{i,j} { mysort() ; }

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
    friend ostream& operator<< (ostream&, const PolyIndx&) ;
    } ;

class PolyTerm : public Term<real,PolyIndx>		// Polynomial term
    {
    public:
    using	Term<real,PolyIndx>::Term ;

    operator	PolyIndx()	   const { return item ; }	// Conversion
    const numb	operator[](uint i) const { return item[i] ; }	// Subscript
    numb&	operator[](uint i)	 { return item[i] ; }	// Subscript

    PolyTerm	operator* (const PolyTerm&) const ;		// Combine
    PolyTerm	operator* (real z) const			// Scale term
		{ PolyTerm ans { *this } ; ans.coeff *= z ; return ans ; }

    int		order()    const { return item.order() ; }	// Term order
    bool	validate() const { return item.validate() ; }	// Validate
    void	print (ostream&, const ObsList&) const ;	// Print
    } ;

struct Polyhash						// PolyIndx hash function
    {
    std::size_t operator() (const PolyIndx&) const noexcept ;
    } ;

class PolyMap : public hash<PolyIndx,real,Polyhash>	// Observable polynomial
    {
    std::reference_wrapper<ObsList> list ;		// Observable list

    public:
    PolyMap (ObsList& l) : list(l) {}			// Constructor

    ObsList& obslist() const { return list ; }		// Underlying ObsList

    void add (PolyIndx t, real d)			// Add term to PolyMap
	{
	auto [iter, isnew] { try_emplace (t,d) } ;
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
	std::erase_if (*this, [](const pair<PolyIndx,real>& p)
	    { return !p.second ; }) ;
	}
    bool allzero () { purge() ; return !size() ; }	// Vanishing Poly?

    friend ostream& operator<< (ostream&, const PolyMap&) ;
    } ;

class ObsPoly : public vector<PolyTerm>			// Polynomial of Obs
    {
    std::reference_wrapper<ObsList> list ;		// Observable list

    public:
    void	push_map (PolyMap&) ;			// PolyMap -> ObsPoly
    void	add (const PolyElem&) ;			// Add term
    void	sort() ;				// Sort terms
    bool	allzero() const ;			// Vanishnig poly?
    ObsPoly&	negate() ; 				// Negate poly
    ObsList&	obslist() const { return list ; }	// Underlying ObsList

    ObsPoly (ObsList& l) : list(l) {}			// Constructor
    ObsPoly (numb indx, ObsList& l)			// Constructor
	    : list(l), vector<PolyTerm>(1, PolyTerm(indx)) {}

    friend ostream& operator<< (ostream&, const ObsPoly&) ;
    } ;

class PolyElem : public Element				// Packed ObsPoly
    {
    public:
    const PolyElem* begin() const { return this + 1 ; }
    const PolyElem* end()   const { return this + this->hdr.len + 1 ; }

    static PolyTerm nextterm (const PolyElem*& ptr)
	{
	PolyTerm term { PolyIndx{}, ptr++->coeff } ;
	for (int i(0) ; i < PSIZ ; ++i)
	    {
	    numb indx { ptr++->index } ;
	    term.item[i] = abs(indx) ;
	    if (indx >= 0) break ;
	    }
	return term ;
	}

    friend ostream& operator<< (ostream&, const PolyElem&) ;
    } ;

static_assert (sizeof(PolyElem) == sizeof(numb)) ;

class PolyRec : public DataRec				// Polynomial data record
    {
    public:
    mutable vector<ulong> offset ;			// Item offests

    using DataRec::DataRec ;

    void add (PolyMap&) ;				// Add PolyTerms
    void add (const ObsPoly&) ;				// Add PolyTerms

    void clear () { DataRec::clear() ; offset.clear() ; }	// Clear data

    const PolyElem& operator() (numb k, numb j=0, numb i=0) const // Indexed Poly
	{						// N.B. slice/col major!
	if (offset.empty()) slice_n_dice () ;
	auto dataptr { begin().ptr } ;
	return dataptr [offset [(i * entry().ncol + j) * entry().nrow + k]] ;
	}

    void slice_n_dice () const				// Initialize offset array
	{
	auto n { entry().items() } ;
	auto m { n } ;

	offset.reserve (n) ;
	for (const auto& poly : *this)
	    {
	    offset.push_back (&poly - begin().ptr) ; --m ;
	    }
	if (m) gripe (format("Corrupted data block: expected {}, got {}",n,n-m)) ;
	}

    struct PolyIter					// Packed PolyElem iterator
	{
	const PolyElem*	ptr ;
	const PolyElem&	operator*()  { return *ptr ; }
	PolyIter& 	operator++() { ptr += ptr->hdr.len + 1 ; return *this ; }
	bool		operator!= (const PolyIter& a) { return ptr != a.ptr ; }

	PolyIter (const PolyElem* p) : ptr(p) {}
	} ;

    const PolyIter begin() const
	{ return PolyIter (static_cast<const PolyElem*>(data())) ; }

    const PolyIter end()   const
	{ return PolyIter (static_cast<const PolyElem*>(data()) + size()) ; }

    } ;

template <size_t N>
class PolyArr : public array<PolyRec,N>			// Array of PolyRec's
    {
    public:
    PolyArr (SysIndex& indx, RecordID id)		// Public constructor
	: PolyArr (std::make_index_sequence<N>{}, indx, id) {}

    private:
    template <size_t... Is> constexpr			// Private constructor
    PolyArr (std::index_sequence<Is...>, SysIndex&, RecordID) ;
    } ;

#endif
