#ifndef ASSESS_H
#define ASSESS_H
#include "Theory.h"

class Site
    {
    public:
    ulong	point ;		// site coordinate
    short	pos ;		// position in SymbStr
    short	nEs ;		// # preceeding E's
    symb	in ;		// prior Symb
    symb	out ;		// following Symb

    bool operator== (const Site& b) const noexcept
	{
	return point == b.point ;
	} ;
    bool operator< (const Site& b) const noexcept
	{
	return point < b.point ;
	} ;
    } ;

class Intersect
    {
    public:
    float	score ;		// factorization score
    short	pos1 ;		// start of subloop1
    short	pos2 ;		// start of subloop2
    short	nEs ;		// # intervening E's

    static constexpr int maxtry { 20 } ;	// max approx tries

    bool operator< (const Intersect& b) const noexcept
	{
	return score < b.score ;
	} ;
    } ;

#endif

