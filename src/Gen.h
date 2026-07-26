#ifndef GEN_H
#define GEN_H
#include "Obs.h"
#include "Op.h"

class Proj ;

class Gen : public OpSum		// Coherence algebra generator
    {
    public:
    short	order { -1 } ;			// Generator order
    OpType	type { OpType::Invalid } ;	// Generator type
    bool	T_odd ;				// Time reversal odd?
    bool	imag ;				// Imaginary coefficient?
    bool	active { true } ;		// Active generator?
    real	coeff { 1 } ;			// Overall coefficient
    ObsPoly	reduction { ObsList::redu } ;	// "inner" commutator

    Gen (OpList& l) : OpSum (l) {}		// Constructor
    Gen (const Op&) ;				// Constructor
    Gen (const Op&,    const Proj&) ;		// Constructor
    Gen (const OpSum&, const Proj&) ;		// Constructor
    Gen (const OpSum&, Element&, real) ;	// Constructor

    Gen&	collect		() ;		// Collect terms
    void 	inner_commute	() ;		// Evaluate reduction
    void 	settype		(Op&) ;		// Set type, order
    void	normalize	(int) ;		// Normalize generator
    bool	allzero		() const ;	// Vanishing Gen?
    ostream&	print (ostream&) const ;	// Short print form

    bool valid()	const { return type != OpType::Invalid ; }
    bool is_Loop()	const { return type == OpType::Loop ; }
    bool is_Eloop()	const { return type == OpType::Eloop ; }
    bool is_Fermion()	const { return type == OpType::Fermion ; }
    bool isgauge()	const { return type == OpType::Loop ||
				       type == OpType::Eloop ; }

    static inline bool autoToddgens { true } ;	// Use commutator E-gens?
    static inline bool gennorm   { false } ;	// Normalize generators?

    static bool	isnew (int, const Gen&) ;	// Dependency test
    static int	project		(Op&) ;		// Project Op onto reps
    static int	project		(OpSum&&) ;	// Project Op sum
    static int	addgen		(OpSum&&) ;	// Add generator
    static void	suspend_group	(uint) ;	// Suspend gen group
    static void	activate_group	(uint) ;	// Activate gen group
    static void	suspend_gen	(uint) ;	// Suspend generator
    static void	activate_gen	(uint) ;	// Activate generator
    static void	activate_gen	() ;		// Activate generators
    static void	normalize	() ;		// Normalize all Gen's
    static void	geninit		(int) ;		// Initialization

    friend ostream& operator<< (ostream&, const Gen&) ;
    } ;

#endif
