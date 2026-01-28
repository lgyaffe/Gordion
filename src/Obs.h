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

class Obs : public SymbStr			// Obs = SymbStr + meta-info
    {
    public: 
    short		corder { -1 } ;	// Creation order
    mutable short	xorder { -1 } ;	// Expectation order
    mutable short	midE   { -1 } ;	// mid-E location
    ObsType		type ;		// Observable type

    explicit Obs (const Op&) ;					// Constructor
    explicit Obs (const string&) ;				// Constructor
    explicit Obs (const string&, ObsType, short, short) ;	// Constructor

    Obs(const SymbStr& s, ObsHdr& hdr)				// Constructor
	:
	SymbStr	(s),
	type	((ObsType) hdr.type),
	corder	(hdr.corder),
	xorder	(hdr.xorder)
	{}

    Obs(const SymbStr& s, ObsType t, short cord, short xord)	// Constructor
	:
	SymbStr (s),
	type	(t),
	corder	(cord),
	xorder	(xord)
	{ if (Obs::check) validate() ; }

    Obs(const_iterator beg, const_iterator end, ObsType t, short cord=-1, short xord=-1)
	:
	SymbStr	(string (beg,end)),
	type	(t),
	corder	(cord),
	xorder	(xord)
	{ if (Obs::check) validate() ; } ;			

    int			canon() ;			// Canonicalize
    int			findstart() noexcept ;		// Rotate to start
    short		middleE()   const ;		// Return mid-E loc
    bool		Esublat()   const ; 		// E sublattice?
    void		validate()  const ;		// Validate
    vector<Site>&	sitelist()  const ;		// Sort site coords

    int		trans	    (const Symm&, int) noexcept;// Symmetry transform
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
    bool is_FermionE()	const { return is_Fermion() && !(size() % 2) ; }
    bool is_FermionO()	const { return is_Fermion() &&  (size() % 2) ; }
    bool imag()		const { return theory.euclid && is_FermionE() ; }

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
/*
    bool is_coord() const { return is_Loop()
				|| is_EEloop()
				|| is_Efermion() && oddlen() ^ !staggered()
				|| is_Fermion() &&
				    (theory.euclid && !staggered() ||
				    !theory.euclid && staggered() ^ oddlen()) ; }
*/
    bool is_coord() const { return is_Loop()
				||  theory.euclid && is_Fermion() && !staggered()
				|| !theory.euclid && (is_EEloop() ||
				    is_Fermion() && staggered() ^ oddlen()) ; }

    bool oddlen() const					// Odd length bilinear?
	{
	return (size() % 2) ^ (is_Efermion() && nostep((*this)[middleE()])) ;
	}

    static inline bool	check { false } ;		// Extra validity tests?
    static ObsType	obstype (const string) ;	// Determine Obs type

    static constexpr const char* type_name[] = 		// ObsType names
	    { "Loop", "Eloop", "Fermion", "EEloop", "Efermion", "Entropy" } ;

    friend ostream& operator<< (ostream&, const Obs&) ;

    static inline const std::map<string,ObsType> obstypes
	{
	{ "Loop",     ObsType::Loop},
	{ "Eloop",    ObsType::Eloop},
	{ "EEloop",   ObsType::EEloop},
	{ "Fermion",  ObsType::Fermion},
	{ "Efermion", ObsType::Efermion},
	{ "Entropy",  ObsType::Entropy}
	} ;
    } ;

struct Obshash					// Obs hash function 
    {
    std::size_t operator()(const SymbStr& s) const
	{
	return std::hash<string>{}(s) ;
	}
    using is_transparent = void ;
    } ;

struct Obs_eq					// Obs equality function 
    {
    bool operator()(const SymbStr& s, const SymbStr& t) const
	{
	return std::equal_to<string>{}(s,t) ;
	}
    using is_transparent = void ;
    } ;

using Obsset = unordered_set<Obs,Obshash,Obs_eq> ;
using Obsmap = hash<Obs,uint,Obshash,Obs_eq> ;

class ObsList: public vector<const Obs*>
    {
    public:
    Obsmap 		map ;			// Hash table
    string		name ;			// List name
    bool		canonicalize ;		// Canonicalize entries?
    bool		classify ;		// Classify entries?
    bool		approx {false} ;	// Approximate exclusions?

    ObsList (const string, bool=false, bool=false) ;

    const Obs&	operator() (uint indx) const	// Return indexed Obs
		    { return *(*this)[indx] ; }

    uint	find (const SymbStr& s) const	// Find Obs, return index
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
    void	purge (uint limit)		// Purge entries
		    {
		    resize (limit) ;
		    std::erase_if (map, [limit](const auto& p)
			{ return p.second >= limit ; }) ;
		    }
    void	insert	(Obsset& set)		// Merge new Obs set
		    {
		    reserve (size() + set.size()) ;
		    map.reserve (size() + set.size()) ;
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

    uint nobs ()	   	 const { return size() ; }
    bool neq  (const ObsList& l) const { return this != &l ; }

    void	obsinit  () ;			// Load basic Obs
    void	do_fermi_init() ;		// Fermion -> Loop map
    uint	store	 (const Obs&) ;		// Store in list
    PolyTerm	catalog  (Obs) ;		// Catalog Obs
    PolyTerm	catalog  (Obs, Obs) ;		// Catalog Obs
    PolyTerm	assess   (Obs&) ;		// Store, approx or discard?
    PolyTerm	is_known (Obs&&) const ;	// Find in list
    PolyTerm	is_known (Obs&&, Obs&&) const ;	// Find in list

    ostream& print(ostream&, uint) const ;	// Print obs
    ostream& print(ostream&) const ;		// Print list

    static ObsList			obs  ;		// Canonicalized Obs
    static ObsList			base ;		// Basic defined Obs
    static ObsList			redu ;		// Gen reductions
    static inline thread_local Obsset	inbox ;		// Obs awaiting insertion
    static inline thread_local bool	freeze {true} ;	// Freeze master list?
    static inline vector<uint2>		fermiinit ;	// Fermion -> Loop map

    static void	retain (const Obs& o)		// Retain for later insertion
	{ inbox.insert (o) ; }
    } ;

class ObsSet : public std::set<uint>				// Obs index set
    {
    public:
    friend ostream& operator<< (ostream& stream, const ObsSet& s)
	{
	for (auto i : s) stream << ObsList::obs(i) << " " ;
	return stream ;
	}
    } ;

class ObsStats : public array<vector<vector<uint>>,nobstype>	// Obs statistics
    {
    public:

    ObsStats (const ObsList&) ;

    uint get(int t, int c, int x) const
	{
	return (t < size()
		&& c >= 0 && c < (*this)[t].size()
		&& x >= 0 && x < (*this)[t][c].size())
		? (*this)[t][c][x] : 0 ;
	}
    uint get(int c, int x) const
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
    uint  maxloop ;
    } ;

#endif
