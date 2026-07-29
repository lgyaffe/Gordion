#ifndef OBS_H
#define OBS_H
#include "Symb.h"
#include "Poly.h"
#include "Assess.h"
#include <climits>
#include <map>

enum class ObsType : char			// Observable categories
    {
    Loop,
    Eloop,
    Fermion,
    EEloop,
    Efermion,
    Entropy,
    _Count_
    } ;
constexpr int nobstype { (int) ObsType::_Count_ } ;

class Op ;
class Symm ;
class ObsList ;

class Obs : public Str				// Obs = Str + meta-info
    {
    public: 
    short		corder { -1 } ;		// Creation order
    mutable short	xorder { -1 } ;		// Expectation order
    mutable short	midE   { -1 } ;		// mid-E location
    ObsType		type ;			// Observable type

    explicit Obs (const Op&) ;					// Constructor
    explicit Obs (const string&) ;				// Constructor
    explicit Obs (const string&, ObsType, short, short) ;	// Constructor

    Obs(const Str& s, const ObsHdr& hdr)			// Constructor
	:
	Str	(s),
	type	((ObsType) hdr.type),
	corder	(hdr.corder),
	xorder	(hdr.xorder)
	{}

    Obs(const Str& s, ObsType t, short cord, short xord)	// Constructor
	:
	Str	(s),
	type	(t),
	corder	(cord),
	xorder	(xord)
	{ if (check) validate() ; }

    Obs(const_iterator beg, const_iterator end, ObsType t, short cord=-1, short xord=-1)
	:
	Str	(beg,end),
	type	(t),
	corder	(cord),
	xorder	(xord)
	{ if (check) validate() ; } ;			

    int		canon() ;				// Canonicalize
    int		findstart() noexcept ;			// Rotate to start
    short	middleE()   const ;			// Return mid-E loc
    bool	Esublat()   const ; 			// E sublattice?
    void	validate()  const ;			// Validate
    SiteVec&	sitelist()  const ;			// Sort site coords
    int		trans (const Symm&, int) noexcept;	// Symmetry transform

    PolyTerm	approximate (const ObsList&) const ;	// Factor/approx
    short	noEbound    (const ObsList&) const ;	// no-E xorder bound
    bool	classify    (ObsList&) ;		// Determine xorder
    short	u1bound     (short) const ;		// U(1) xorder bound
    bool	reduce      (short, ObsList&) ;		// Commute w. plaq/hop
    bool	factorize   (short, const ObsList&) const ;	// Factorize

    short order()	const { return corder + xorder ; }
    bool  known_xord()	const { return xorder >= 0 ; }
    bool  bilinear()	const { return type == ObsType::Fermion ||
				       type == ObsType::Efermion ; }
    bool  EorEEloop()	const { return type == ObsType::Eloop ||
				       type == ObsType::EEloop ; }
    bool  midEtype()	const { return type == ObsType::EEloop ||
				       type == ObsType::Efermion ; }
    bool staggered()	const { return bilinear() && isstag(front())
						   ^ isstag(back()) ; }

    bool is_Loop()	const { return type == ObsType::Loop ; }
    bool is_Eloop()	const { return type == ObsType::Eloop ; }
    bool is_EEloop()	const { return type == ObsType::EEloop ; }
    bool is_Fermion()	const { return type == ObsType::Fermion ; }
    bool is_Efermion()	const { return type == ObsType::Efermion ; }
    bool is_Entropy()	const { return type == ObsType::Entropy ; }
    bool is_FermionO()	const { return is_Fermion() &&  (size() % 2) ; }
    bool imag()		const { return theory.euclid && !is_FermionO() ; }

    bool is_gauge() const { return is_Loop()
				|| is_Eloop()
				|| is_EEloop() && !isEE_F(front())
				|| front() == EntrG ; }
    bool has_Es()   const { return is_Eloop()
				|| is_EEloop()
				|| is_Efermion() ; }
    bool is_fermi() const { return is_Fermion()
				|| is_Efermion()
				|| isEE_F(front())
				|| front() == EntrF ; }
    bool is_coord() const { return theory.euclid
				    ? is_Loop() || is_Fermion() && !staggered()
				    : is_Loop() || is_EEloop() ||
				      is_Fermion()  && staggered() ^  oddlen() ||
				      is_Efermion() && staggered() ^ !oddlen() ; }

    bool oddlen() const					// Odd length bilinear?
	{
	return (size() % 2) ^ (is_Efermion() && nostep((*this)[middleE()])) ;
	}

    static inline bool	check { false } ;		// Extra validity tests?
    static ObsType	obstype (const string) ;	// Determine Obs type

    static constexpr const char* type_name[] = 		// ObsType names
	    { "Loop", "Eloop", "Fermion", "EEloop", "Efermion", "Entropy" } ;

    static inline const std::map<string,ObsType> obstypes
	{
	{ "Loop",     ObsType::Loop},
	{ "Eloop",    ObsType::Eloop},
	{ "EEloop",   ObsType::EEloop},
	{ "Fermion",  ObsType::Fermion},
	{ "Efermion", ObsType::Efermion},
	{ "Entropy",  ObsType::Entropy}
	} ;

    friend ostream& operator<< (ostream&, const Obs&) ;
    } ;

using Obsset = unordered_set<Obs,Strhash,Str_eq> ;
using Obsmap = hash<Obs,numb,Strhash,Str_eq> ;

class ObsList: vector<const Obs*>
    {
    Obsmap 		map ;			// Hash table
    public:
    string		name ;			// List name
    bool		canonicalize ;		// Canonicalize entries?
    bool		classify ;		// Classify entries?
    bool		approx {false} ;	// Approximate exclusions?

    using vector::const_iterator ;	// Expose base const_iterator
    const_iterator begin() const { return vector::cbegin() ; }
    const_iterator end()   const { return vector::cend()   ; }

    ObsList (const string, bool=false, bool=false) ;

    const Obs&	operator() (numb indx) const	// Return indexed Obs
		    { return *(*this)[indx] ; }

    auto	size() const { return vector::size() ; } // List size
    auto	shrink() { return shrink_to_fit() ; }	// Shrink

    numb	find (const Str& s) const	// Find Obs, return index
		    {
		    auto p1 { map.find(s) } ;
		    if (p1 != map.end()) return p1->second ;
		    if (this == &ObsList::obs && inbox.size())
			{
			auto p2 { inbox.find(s) } ;
			if (p2 != inbox.end()) return UINT_MAX-1 ;
			}
		    return UINT_MAX ;
		    }
    void	reserve (int len)		// Reserve space
		    {
		    vector::reserve (len) ;
		    map.reserve (len) ;
		    }
    void	insert	(Obsset& set)		// Merge new Obs set
		    {
		    reserve (size() + set.size()) ;
		    for (const auto& obs : set) store (obs) ;
		    }
    bool	frozen () const			// Frozen master list?
		    {
		    return this == &ObsList::obs && ObsList::freeze ;
		    }
    bool	freezeif () const		// Freeze if master list
		    {
		    if (this == &ObsList::obs)
			{
			bool prev { ObsList::freeze } ;
			ObsList::freeze = true ;
			return prev ;
			}
		    return false ;
		    }
    void	refreezeif (bool prev) const	// Reset freeze if master list
		    {
		    if (this == &ObsList::obs)
			{
			ObsList::freeze = prev ;
			}
		    }
    bool	 neq  (const ObsList& l) const { return this != &l ; }

    int		do_fermiinit () ;		// Load Fermion -> Loop map
    void	obsinit	(int) ;			// Load basic Obs
    ostream&	print	(ostream&, numb) const ;// Print obs
    ostream&	print	(ostream&) const ;	// Print list

    numb	store	 (const Obs&) ;		// Store in list
    void	purge	 (numb)	;		// Purge entries
    void	clear	 () ;			// Clear list
    void	rehash	 () const ; 		// Recalculate hashes
    void	hasher	 (ulong&,const Obs&) const ; // List hasher
    PolyTerm	catalog  (Obs) ;		// Catalog Obs
    PolyTerm	catalog  (Obs, Obs) ;		// Catalog Obs
    PolyTerm	is_known (Obs&&) const ;	// Find in list
    PolyTerm	is_known (Obs&&, Obs&&) const ;	// Find in list
    PolyTerm	assess   (Obs&) ;		// Store, approx or discard?

    static ObsList			obs  ;		// Canonicalized Obs
    static ObsList			base ;		// Basic defined Obs
    static ObsList			redu ;		// Gen reductions
    static inline vector<numb2>		fermiinit ;	// Fermion -> Loop map
    static inline thread_local Obsset	inbox ;		// Obs awaiting insertion
    static inline thread_local bool	freeze {true} ;	// Freeze master list?
    static inline bool			swapped {false};// Swapped to disk?

    static void	ondisk () ;			// Leave ObsList::obs on disk
    static void	retain (const Obs& o)		// Retain for later insertion
	{ inbox.insert (o) ; }
    } ;

class ObsSubset : public std::map<numb,Obs>		// Obs subset
    {
    public:
    friend ostream& operator<< (ostream& stream, const ObsSubset& s)
	{
	for (auto& [i,o] : s) stream << o << " " ;
	return stream ;
	}
    } ;

class ObsStats : public array<vector<vector<ulong>>,nobstype>	// Obs statistics
    {
    public:

    ObsStats (const ObsList&) ;

    ulong get(int t, int c, int x) const
	{
	return (t < size()
		&& c >= 0 && c < (*this)[t].size()
		&& x >= 0 && x < (*this)[t][c].size())
		? (*this)[t][c][x] : 0 ;
	}
    ulong get(int c, int x) const
	{
	int sum (0) ;
	for (int t(0) ; t < size() ; ++t)
	    {
	    if (c >= 0 && c < (*this)[t].size() &&
		x >= 0 && x < (*this)[t][c].size())
		sum += (*this)[t][c][x] ;
	    }
	return sum ;
	}
    short maxc = 0 ;
    short maxx = 0 ;
    float avglen ;
    numb  maxloop ;
    } ;

#endif
