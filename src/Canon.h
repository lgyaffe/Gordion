#ifndef CANON_H
#define CANON_H
#include "Global.h"
#include "Index.h"
#include "Obs.h"

static constexpr int loopchunksize = theory.dim == 4 ? 6 :	// Loop chuck size
				     theory.dim == 3 ? 7 : 9 ;

static constexpr int specchunksize = !theory.euclid		// Spec chunk size
				   ? theory.dim == 2 ? 5 : 4
				   : theory.dim == 4 ? 6 :
				     theory.dim == 3 ? 7 : 8 ;

static constexpr int loopbits = 3 ;				// Loop score bits
static constexpr int specbits = theory.euclid ? 3 : 6 ;		// Spec score bits

static constexpr ulong looptblsize (int len = loopchunksize)	// Loop table size
    {
    return ipow(2 * theory.dim, len) ;
    }

static constexpr ulong spectblsize (int len = specchunksize)	// Special table size
    {
    return theory.euclid ?
	   2 * ipow (2 * theory.dim, len) * !!theory.nf :
	   2 * ipow (2 * theory.dim, len) * (4 + 2 * !!theory.nf) * len ;
    }

static_assert (ipow(2,loopbits*loopchunksize) < UINT_MAX, "loopchunksize too big") ;
static_assert (ipow(2,specbits*specchunksize) < UINT_MAX, "specchunksize too big") ;

class CanonCache : public hash<string,int>		// Short Obs cache
    {
    public:
    bool freeze { false } ;

    void purge (uint indx)		// purge entries >= indx
	    { std::erase_if (*this, [indx](const auto& p)
		{ return abs(p.second) >= indx ; }) ; }

    void load (const Obsset&) ;		// load short inbox Obs into cache
    void reload () ;			// reload cache

    friend ostream& operator<< (ostream&, const CanonCache&) ;
    } ;

using SymmSet = vector<uint> ;	 			// Trial symmetry subset

class Node						// Canonicalization chunk data
    {
    public:
    uint	score = UINT_MAX ;	// chunk score
    int		indx = -1 ;		// indx of symmetry subset to try
    } ;

class Canon						// Canonicalization tables
    {
    friend Obs ;

    static uint speccode(const Obs&, int, bool) ;
    static bool loopchunk(uint, SymbStr&) ;
    static bool specchunk(uint, SymbStr&) ;

    public:
    static inline vector<Node>		looptable { looptblsize() } ;
    static inline vector<Node>		spectable { spectblsize() } ;
    static inline Index<SymmSet>	symmset ;
    static inline CanonCache		cache ;
    static inline struct Stats
	{
	uint	symmsubsets  ;
	uint	looptblnodes ;
	uint	spectblnodes ;
	}				stats ;

    static void looptblinit() ;
    static void spectblinit() ;
    } ;

#endif
