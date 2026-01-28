#ifndef THEORY_H
#define THEORY_H
#include "Gordion.h"
#include "Gripe.h"

enum LATTTYPE { Hamilton = 0, Euclid = 1 } ;	// Basic lattice theory class

class Coord					// Lattice coordinate vector
    {
    public:
    union { short comp[4] ; ulong point ; } ;

    void modlen (int dir, Coord box)		// put single component into box
	{
	if (box.comp[dir])
	    {
	    while (comp[dir] < 0)		comp[dir] += box.comp[dir] ;
	    while (comp[dir] >= box.comp[dir])	comp[dir] -= box.comp[dir] ;
	    }
	}
    void modulo (Coord box)			// put all comp's into box
	{
	for (int i(0) ; i < 4 ; ++i) modlen (i,box) ;
	}
    bool isclosed (Coord box)			// vanishing displacement?
	{
	if (box.point) modulo (box) ;
	return !point ;
	}
    } ;

class Theory					// Basic theory info
    {
    public:
    const char8	 	name ;			// theory name
    const char		dim ;			// lattice dimension
    const char		nf  ;			// number of fermions
    const int		nrep ;			// number of irreps
    const bool		euclid { false } ;	// Hamiltonian or Euclidean?
    const Coord		box { 0,0,0,0,} ;	// lattice size: 0:=infinity

    constexpr bool unbounded() const		// unbounded lattice?
	{ return !box.comp[0] && !box.comp[1] &&
		 !box.comp[2] && !box.comp[3] ; }

    constexpr bool oddboxsize() const		// odd length lattice size?
	{ return box.comp[0] % 2 || box.comp[1] % 2 ||
		 box.comp[2] % 2 || box.comp[3] % 2 ; }

    string box_size() const			// pretty print box dimensions
	{
	string ans ;
	for (int d(0) ; d < dim ; ++d)
	    {
	    uint len ( box.comp[d] ) ;
	    ans += len ? format ("{} ", len) : "\u221E " ;
	    }
	return ans ;
	} ;

    char8 parent() const			// YM parent theory
	{
	string	tmp { name.begin(), name.end() } ;
	if (tmp.find ("QCD",0) == 0)
	    {
	    char8 YMthy ;
	    tmp.replace (0,3,"YM") ;
	    std::copy (tmp.begin(), tmp.end(), YMthy.begin()) ;
	    return YMthy ;
	    }
	else return name ;
	}

    constexpr bool isR1() const { return dim == 1 && unbounded() ; }
    constexpr bool isR2() const { return dim == 2 && unbounded() ; }
    constexpr bool isR3() const { return dim == 3 && unbounded() ; }
    constexpr bool isR4() const { return dim == 4 && unbounded() ; }

    constexpr bool isS1()    const { return dim == 1 &&  box.comp[0] ; }
    constexpr bool isR1xS1() const { return dim == 3 && !box.comp[0] &&  box.comp[1] ; }
    constexpr bool isR2xS1() const { return dim == 3 && !box.comp[0] && !box.comp[1]
						     &&  box.comp[2] ; }
    constexpr bool isR3xS1() const { return dim == 4 && !box.comp[0] && !box.comp[1]
						     && !box.comp[2] &&  box.comp[3] ; }

    static void theoryinit() ;
    static void theorydefn() ;
    } ;

#if   defined (YM1h)
    constexpr Theory theory { "YM1h",  1, 0, 4,  Hamilton, { 1 } } ;
#elif defined (YM2h)
    constexpr Theory theory { "YM2h",  2, 0, 12, Hamilton } ;
#elif defined (YM3h)
    constexpr Theory theory { "YM3h",  3, 0, 40, Hamilton } ;
#elif defined (YM1e)
    constexpr Theory theory { "YM1e",  1, 0, 4,  Euclid, { 1 } } ;
#elif defined (YM2e)
    constexpr Theory theory { "YM2e",  2, 0, 12, Euclid } ;
#elif defined (YM3e)
    constexpr Theory theory { "YM3e",  3, 0, 40, Euclid } ;
#elif defined (YM4e)
    constexpr Theory theory { "YM4e",  4, 0, 2,  Euclid } ;
#elif defined (QCD1h)
    constexpr Theory theory { "QCD1h", 1, 1, 4,  Hamilton } ;
#elif defined (QCD2h)
    constexpr Theory theory { "QCD2h", 2, 1, 12, Hamilton } ;
#elif defined (QCD3h)
    constexpr Theory theory { "QCD3h", 3, 1, 40, Hamilton } ;
#elif defined (QCD4h)
    constexpr Theory theory { "QCD4h", 4, 1,  2, Hamilton } ;
#elif defined (QCD1e)
    constexpr Theory theory { "QCD1e", 1, 1,  4, Euclid } ;
#elif defined (QCD2e)
    constexpr Theory theory { "QCD2e", 2, 1, 12, Euclid } ;
#elif defined (QCD3e)
    constexpr Theory theory { "QCD3e", 3, 1, 40, Euclid } ;
#elif defined (QCD4e)
    constexpr Theory theory { "QCD4e", 4, 1, 2,  Euclid } ;
#elif defined (YM3hp)
    constexpr Theory theory { "YM3hp", 3, 0, 24, Hamilton, { 0, 0, 1 } } ;
#elif defined (YM3ep)
    constexpr Theory theory { "YM3ep", 3, 0, 24, Euclid,   { 0, 0, 1 } } ;
#else
#error "No known theory selected!"
#endif

static_assert (theory.dim >= 1 && theory.dim <= 4, "Invalid theory dimension") ;
static_assert (!theory.nf || !theory.oddboxsize(), "Need even lattice period(s)") ;
static_assert (theory.nf <= 2, "Too many fermion flavors") ;

static constexpr int PSIZ = theory.euclid ? 2 : 4 ;	// Polyterm order
static constexpr int NREP = theory.nrep ;		// Sum of irrep dimensions

#endif
