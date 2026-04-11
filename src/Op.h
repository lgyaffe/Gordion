#ifndef OP_H
#define OP_H
#include "Symb.h"
#include "Obs.h"
#include "Term.h"
#include "Index.h"
#include "Global.h"

using OpTerm = Term<doub,uint> ;			// Term with list index

class OpSum : public vector<OpTerm>			// Op linear combination
    {
    public:
    using	vector<OpTerm>::vector ;
    int		collect (bool = false) ;		// collect terms

    friend ostream& operator<< (ostream&, const OpSum&) ;
    } ;

enum class OpType : char				// Operator types
    {
    Loop,
    Eloop,
    Fermion,
    Invalid
    } ;

class Op : public SymbStr				// Op = generator term
    {
    public:
    OpType		type { OpType::Invalid } ;	// Observable type
    mutable bool	primary { false } ;		// Primary Op?
    short		order { -1 } ;			// Grading order

    Op () {}						// Default constructor
    explicit Op (const string&, short) ;		// Constructor
    explicit Op (const string&, OpType, short) ;	// Constructor

    Op (const SymbStr& s, OpHdr& hdr)			// Constructor
	:
	SymbStr	(s),
	type	((OpType) hdr.type),
	order	(hdr.order),
	primary	(hdr.prim)
	{ if (Op::check) validate() ; }

    Op (const SymbStr& s, OpType t, short o)		// Constructor
	:
	SymbStr	(s),
	type	(t),
	order	(o)
	{ if (type == OpType::Loop) findstart() ;
	  if (Op::check) validate() ; }

    Op (const Obs& a) 					// Convert Obs -> Op
	:
	SymbStr	(a),
	order	(a.corder),
	type	(a.is_Fermion() ? OpType::Fermion :
		 a.is_Eloop() ? OpType::Eloop :
		 a.is_Loop() ? OpType::Loop : OpType::Invalid)
	{}

    void findstart() ;
    bool oddlen()	const { return size() % 2 ; }
    bool is_Loop()	const { return type == OpType::Loop ; }
    bool is_Eloop()	const { return type == OpType::Eloop ; }
    bool is_Fermion()	const { return type == OpType::Fermion ; }
    bool is_FermionE()	const { return is_Fermion() && !oddlen() ; }
    bool is_FermionO()	const { return is_Fermion() &&  oddlen() ; }
    bool staggered()	const { return is_Fermion() &&
				isstag(front()) ^ isstag(back()) ; }
    bool is_coord()	const { return is_Loop() || is_Fermion() &&
				(theory.euclid && !staggered() ||
				!theory.euclid && staggered() ^ oddlen()) ; }

    static void		setprimary () ;			// Identify primary Ops
    static OpSum 	loop_dt (Op) ;			// Return [EE,loop]/2
    static OpSum 	loop_dt (OpSum) ;		// Return [EE,loops]/2
    static OpSum 	flipT (OpSum) ;			// Flip fermion staggering
    static OpType	optype (const string) ;		// Determine type
    static uint		store (const Op&) ;		// Store in list
    static void		purge (uint limit)		// Purge entries
	    {
	    auto& nopG { global.nopG() } ;
	    auto k { std::erase_if (list.map, [limit,nopG](const auto& p)
		{ return p.second >= limit && p.second >= nopG; }) } ;
	    global.nopF() -= k ;
	    auto l { std::erase_if (list.map, [limit,nopG](const auto& p)
		{ return p.second >= limit && p.second < nopG; }) } ;
	    global.nopG() -= l ;
	    list.resize (limit) ;
	    }

    static inline bool	check { false } ;		// Do validity tests?
    static Index<Op>	list ;				// List of defined Ops

    static ostream& print(ostream&, uint) ;		// Print indexed Op
    static ostream& print(ostream&) ;			// Print Op list
    friend ostream& operator<< (ostream&, const Op&) ;	// Print Op

    private:
    void		validate() ;
    static OpSum 	loop_dt (OpTerm,OpSum&) ;	// Return [EE,loop]/2
    } ;

#endif
