#ifndef LIST_H
#define LIST_H
#include "Hash.h"
#include "Obs.h"
#include "Op.h"
#include "Poly.h"
#include "Str.h"

template <typename T>
class Index : public vector<T>		// vector<T> plus string keyed index
    {
    public:
    hash<string,uint>  map ;			// string key -> index map

    numb store (const string& key, const T& t)	// store if new, return indx
	{
	auto [iter, isnew] { map.try_emplace (key, (*this).size()) } ;
	if (isnew) vector<T>::push_back (t) ;
	return iter->second ;
	}
    void clear ()				// Clear index'd list
	{
	map.clear () ;
	vector<T>::clear () ;
	}
    } ;

template <typename T>
class List : protected vector<const T*>		// Indexed list of T*'s
    {
    protected:
    using ListMap = hash<T,numb,Strhash,Str_eq> ;
    ListMap	map ;					// Hash table

    public:
    using const_iterator = typename std::vector<const T*>::const_iterator ;

    const_iterator	cbegin() const noexcept
			    { return vector<const T*>::cbegin() ; }
    const_iterator	cend() const noexcept
			    { return vector<const T*>::cend() ; }

    auto	size() const				// List length
		    { return vector<const T*>::size() ; }
    const T&	operator() (numb indx) const		// Return indexed T
		    { return *(*this)[indx] ; }

    numb	find (const Str& s) const		// Find T, return index
		    {
		    auto p { map.find(s) } ;
		    if (p != map.end()) return p->second ;
		    return MAXNUM ;
		    }
    void	reserve (int len)			// Reserve space
		    {
		    vector<const T*>::reserve (len) ;
		    map.reserve (len) ;
		    }
    } ;

class OpList : public Index<Op>			// Indexed list of Op's
    {
    public:
    numb	store () ;			// Store Op in list
    void	setprimary () ;			// Identify primary Op's

    numb store (const Op& op)			// Store Op in list
	{
	auto	len  ( size() ) ;
	numb	indx { Index::store (op,op) } ;
	if (indx >= size()) fatal ("OpList::store: bad store! ") ;
	return indx ;
	}
    void purge (numb limit)			// Purge entries
	{
	std::erase_if (map, [&](const auto& p) { return p.second >= limit ; }) ;
	resize (limit) ;
	}
    static void opinit (uint) ;			// List initialization

    ostream& print (ostream&, numb) const ;	// Print indexed Op
    ostream& print (ostream&) const ;		// Print Op list
    } ;

class ObsList: public List<Obs>			// Indexed list of Obs's
    {
    public:
    string		name ;			// List name
    bool		canonicalize ;		// Canonicalize entries?
    bool		classify ;		// Classify entries?
    bool		approx { false } ;	// Approximate exclusions?
    numb		nobsG = 0 ;		// # gauge Obs
    numb		nobsF = 0 ;		// # fermi Obs
    ulong		hashG = 0 ;		// gauge Obs hash
    ulong		hashF = 0 ;		// fermi Obs hash

    ObsList (const string, bool=false, bool=false) ;	// Constructor

    using vector::const_iterator ;		// Base const_iterator
    const_iterator begin() const { return vector::cbegin() ; }
    const_iterator end()   const { return vector::cend()   ; }

    void	obsinit	(uint) ;		// Load basic Obs
    ostream&	print	(ostream&, numb) const ;// Print Obs
    ostream&	print	(ostream&) const ;	// Print list
    bool	neq (const ObsList& l) const { return this != &l ; }

    numb	store	 (const Obs&) ;		// Store in list
    void	purge	 (numb) ;		// Purge entries
    void	clear	 () ;			// Clear list
    void	rehash	 () ;			// Recalculate hashes
    void	hasher	 (ulong&,const Obs&) const ; // List hasher
    PolyTerm	catalog  (Obs) ;		// Catalog Obs
    PolyTerm	catalog  (Obs, Obs) ;		// Catalog Obs
    PolyTerm	is_known (Obs&&) const ;	// Find in list
    PolyTerm	is_known (Obs&&, Obs&&) const ;	// Find in list

    numb&  nobs	(int stage) { return stage ? nobsF : nobsG ; }
    ulong& hash	(int stage) { return stage ? hashF : hashG ; }

    PolyTerm	is_known (Obs&) const ;		// Find in list
    PolyTerm	assess   (Obs&) ;		// Store, approx or discard?
    numb	find (const Str&) const ;	// Find Obs, return index
    void	ondisk	 () ;			// Leave ObsList on disk
    bool	frozen   () const ;		// Frozen master list?
    bool	freezeif () const ;		// Freeze if master list
    void	refreezeif (bool) const ;	// Reset freeze if master list
    int		do_fermiinit () ;		// Load Fermion -> Loop map

    void	insert	(ObsSet& set)		// Merge new Obs set
		    {
		    reserve (size() + set.size()) ;
		    for (const auto& obs : set) store (obs) ;
		    }
    auto	shrink() { return shrink_to_fit() ; } // Shrink

    static void	retain (const Obs&) ;		// Retain for later insertion
    } ;

#endif
