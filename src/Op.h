#ifndef OP_H
#define OP_H
#include "Index.h"
#include "Obs.h"

using OpTerm = Term<real,numb> ;			// Term with list index
class OpList ;
class Obs ;

enum class OpType : char				// Operator types
    {
    Loop,
    Eloop,
    Fermion,
    Invalid
    } ;

class Op : public Str					// Op = generator term
    {
    public:
    OpType	type { OpType::Invalid } ;	// Operator type
    bool	primary { false } ;		// Primary Op?
    short	order   { -1 } ;		// Grading order

    Op () {}						// Default constructor
    explicit Op (const string&, short) ;		// Constructor
    explicit Op (const string&, OpType, short) ;	// Constructor

    Op (const Str& s, const OpHdr& hdr)			// Constructor
	:
	Str	(s),
	type	((OpType) hdr.type),
	order	(hdr.order)
	{ if (check) validate() ; }

    Op (const Str& s, OpType t, short o)		// Constructor
	:
	Str	(s),
	type	(t),
	order	(o)
	{ if (type == OpType::Loop) findstart() ;
	  if (check) validate() ; }

    Op (const Obs& a) 					// Convert Obs -> Op
	:
	Str	(a),
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
    bool is_FermionO()	const { return is_Fermion() &&  oddlen() ; }
    bool staggered()	const { return is_Fermion() &&
				isstag(front()) ^ isstag(back()) ; }
    bool is_coord()	const { return is_Loop() || is_Fermion() &&
				(theory.euclid && !staggered() ||
				!theory.euclid && staggered() ^ oddlen()) ; }

    static OpType	optype (const string) ;		// Determine type
    static inline bool	check { false } ;		// Do validity tests?

    friend ostream& operator<< (ostream&, const Op&) ;	// Print Op

    private:
    void		validate() ;
    } ;

class OpList : public Index<Op>
    {
    public:
    numb	store () ;			// Store Op in list
    void	setprimary () ;			// Identify primary Op's

    numb store (const Op& op)			// Store Op in list
	{
	auto	len  ( size() ) ;
	uint	indx { Index::store (op,op) } ;
	if (indx >= size()) fatal ("OpList::store: bad store! ") ;
	return indx ;
	}
    void purge (numb limit)			// Purge entries
	{
	std::erase_if (map, [&](const auto& p) { return p.second >= limit ; }) ;
	resize (limit) ;
	}

    ostream& print (ostream&, numb) const ;	// Print indexed Op
    ostream& print (ostream&) const ;		// Print Op list
    } ;

class OpSum : public vector<OpTerm>			// Op linear combination
    {
    std::reference_wrapper<OpList> list ;		// Operator list

    public:
    OpSum (OpList&) ;					// Constructor
    OpSum (OpTerm*, OpTerm*, OpList&) ;			// Constructor

    OpList&		oplist() const { return list ;} // Underlying OpList
    int			collect (bool = false) ;	// Collect terms
    OpSum 		flipT   () const ;		// Flip fermion staggering
    OpSum 		loop_dt () ;			// Return [EE,loops]/2
    static OpSum 	loop_dt (Op, OpList&) ;		// Return [EE,loops]/2
    static OpSum 	loop_dt (OpTerm, OpSum&) ;	// Return [EE,loops]/2

    friend ostream& operator<< (ostream&, const OpSum&) ;
    } ;

#endif
