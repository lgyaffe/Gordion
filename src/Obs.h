#ifndef OBS_H
#define OBS_H
#include "Assess.h"
#include "Poly.h"
#include "Str.h"
#include <map>

enum class ObsType : char			// Observable categories
    {
    Loop,
    Eloop,
    Fermion,
    EEloop,
    efermion,
    Entropy,
    _Count_
    } ;
inline constexpr int nobstype { (int) ObsType::_Count_ } ;

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

    explicit Obs (const Op&) ;				// Constructor
    explicit Obs (const string&) ;			// Constructor
    explicit Obs (const string&, ObsType, short,short);	// Constructor

    Obs(const Str& s, const ObsHdr& hdr)		// Constructor
	:
	Str	(s),
	type	((ObsType) hdr.type),
	corder	(hdr.corder),
	xorder	(hdr.xorder)
	{}

    Obs(const Str& s, ObsType t, short cord, short xord) // Constructor
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

    int		canon	  () ;				// Canonicalize
    int		findstart () noexcept ;			// Rotate to start
    short	middleE	  () const ;			// Return mid-E loc
    bool	Esublat	  () const ; 			// E sublattice?
    void	validate  () const ;			// Validate
    SiteVec&	sitelist  () const ;			// Sort site coords
    int		trans (const Symm&, int) noexcept ;	// Symmetry transform

    PolyTerm	approximate (const ObsList&) const ;	// Factor/approx
    short	noEbound    (const ObsList&) const ;	// no-E xorder bound
    bool	classify    (ObsList&) ;		// Determine xorder
    short	u1bound     (short) const ;		// U(1) xorder bound
    bool	reduce      (short, ObsList&) ;		// Commute w. plaq/hop
    bool	factorize   (short, const ObsList&) const ;	// Factorize

    short order()	const { return corder + xorder ; }
    bool  known_xord()	const { return xorder >= 0 ; }
    bool  bilinear()	const { return type == ObsType::Fermion ||
				       type == ObsType::efermion ; }
    bool  EorEEloop()	const { return type == ObsType::Eloop ||
				       type == ObsType::EEloop ; }
    bool  midEtype()	const { return type == ObsType::EEloop ||
				       type == ObsType::efermion ; }
    bool staggered()	const { return bilinear() && isstag(front())
						   ^ isstag(back()) ; }

    bool is_Loop()	const { return type == ObsType::Loop ; }
    bool is_Eloop()	const { return type == ObsType::Eloop ; }
    bool is_EEloop()	const { return type == ObsType::EEloop ; }
    bool is_Fermion()	const { return type == ObsType::Fermion ; }
    bool is_efermion()	const { return type == ObsType::efermion ; }
    bool is_Entropy()	const { return type == ObsType::Entropy ; }
    bool is_FermionO()	const { return is_Fermion() &&  (size() % 2) ; }
    bool is_FermionE()	const { return is_Fermion() && !(size() % 2) ; }
    bool imag()		const { return theory.euclid && is_FermionE() ; }

    bool is_gauge() const { return is_Loop()
				|| is_Eloop()
				|| is_EEloop() && !is_ee(front())
				|| front() == EntrG ; }
    bool has_Es()   const { return is_Eloop()
				|| is_EEloop()
				|| is_efermion() ; }
    bool is_fermi() const { return is_Fermion()
				|| is_efermion()
				|| is_ee(front())
				|| front() == EntrF ; }
    bool is_coord() const { return theory.euclid ?
				is_Loop() || is_Fermion() && !staggered() :
				is_Loop() || is_EEloop() ||
				is_Fermion()  && staggered() ^  oddlen() ||
				is_efermion() && staggered() ^ !oddlen() ; }

    bool oddlen() const				// Odd length bilinear?
	{
	return (size() % 2) ^ (is_efermion() && nostep((*this)[middleE()])) ;
	}
    static ObsType	obstype (const string) ; // Determine Obs type

    static inline bool	check { false } ;	 // Extra validity tests?

    static constexpr const char* type_name[] = 	 // ObsType names
	    { "Loop", "Eloop", "Fermion", "EEloop", "efermion", "Entropy" } ;

    static inline const std::map<string,ObsType> obstypes
	{
	{ "Loop",     ObsType::Loop     },
	{ "Eloop",    ObsType::Eloop    },
	{ "EEloop",   ObsType::EEloop   },
	{ "Fermion",  ObsType::Fermion  },
	{ "efermion", ObsType::efermion },
	{ "Entropy",  ObsType::Entropy  }
	} ;

    friend ostream& operator<< (ostream&, const Obs&) ;
    } ;

using ObsSet = unordered_set<Obs,Strhash,Str_eq> ;

class ObsSubset : public std::map<Obs,numb>		// Obs subset
    {
    public:
    friend ostream& operator<< (ostream& stream, const ObsSubset& s)
	{
	for (const auto& [o,i] : s) stream << o << " " ;
	return stream ;
	}
    } ;

class ObsStats : public array<vector<vector<ulong>>,nobstype>	// Obs statistics
    {
    public:
    ObsStats (const ObsList&) ;

    numb get(int t, int c, int x) const		// given type, corder, xorder
	{
	return (t < size()
		&& c >= 0 && c < (*this)[t].size()
		&& x >= 0 && x < (*this)[t][c].size())
		? (*this)[t][c][x] : 0 ;
	}
    ulong get(int c, int x) const		// given corder, xorder
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
    ulong get(int t) const			// given type
	{
	if (t < size())
	    {
	    const auto& column { (*this)[t] } ;
	    int sum (0) ;
	    for (int c(0) ; c < column.size() ; ++c)
		for (int x(0) ; x < column[c].size() ; ++x)
		    sum += column[c][x] ;
	    return sum ;
	    }
	return 0 ;
	}
    short maxc = 0 ;
    short maxx = 0 ;
    float avglen ;
    numb  maxloop ;
    } ;

#endif
