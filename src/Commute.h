#ifndef COMMUTE_H
#define COMMUTE_H
#include "Poly.h"
#include "Gen.h"

namespace Commute
    {
    void commute_poly (const Gen&, const ObsPoly&, PolyMap&) ;		 // [Gen,Obs]
    void commute_poly (const Gen&, const PolyMap&, PolyMap&) ;		 // [Gen,Obs]
    void commute_term (const Gen&, const PolyTerm&, ObsList&, PolyMap&); // [Gen,Obs]

    Gen& commute_gen (const Gen&, const Gen&,     Gen&) ;		 // [Gen,Gen]
    void op_commute  (doub, const Op&, const Op&, Gen&) ;		 // [Op,Op]

    void do_commute  (const Op&, const Obs&, PolyTerm, ObsList&, PolyMap&) ; // [Op,Obs]
    void do_commuteA (const Op&, const Obs&, PolyTerm, ObsList&, PolyMap&) ; // [Op,Obs]
    void do_commuteB (const Op&, const Obs&, PolyTerm, ObsList&, PolyMap&) ; // [Op,Obs]
    void do_commuteC (const Op&, const Obs&, PolyTerm, ObsList&, PolyMap&) ; // [Op,Obs]
    void do_commuteD (const Op&, const Obs&, PolyTerm, ObsList&, PolyMap&) ; // [Op,Obs]
    void do_commuteE (const Op&, const Obs&, PolyTerm, ObsList&, PolyMap&) ; // [Op,Obs]

    void do_inner (const Obs&, const PolyTerm, ObsList&, PolyMap&) ;     // tr([E,loop])
    void do_split (doub, const ObsPoly&, const Obs&, ObsList&, PolyMap&) ;

    inline short cordsum (const Op& a, const Obs& b)		// Combine corders
	{
	return b.corder >= 0 ? b.corder + a.order : -1 ;
	}

    inline ObsType commtype (const Op& a, const Obs& b)		// [Op,Obs] type, only for
	{							// EEloop or Efermion
	if      (a.is_Loop()    && b.is_EEloop())	return ObsType::Eloop ;
	else if (a.is_Loop()    && b.is_Efermion())	return ObsType::Fermion ;
	else if (a.is_Fermion() && b.is_EEloop())	return ObsType::Efermion ;
	else						return b.type ;
	}
    }

#endif
