#ifndef SYMB_H
#define SYMB_H
#include "Gordion.h"
#include "Theory.h"

class SymbStr : public string				// Symbol string
    {
    public:
    SymbStr	()		   : string() {}
    SymbStr	(const string& s)  : string(s) {}
    SymbStr	(const symb c)     : string(1,c) {}
    SymbStr	(const symb* p, const symb* q) : string(p,q) {}
    SymbStr	(const_iterator p,const_iterator q) : string(p,q) {}

    string	print () const ;
    string	print (const_iterator, const_iterator) const ;
    int		join (symb) ;	
    int		join (const_iterator, const_iterator) ;	
    int		joinends (SymbStr::iterator, SymbStr::iterator) noexcept ;
    void	excise   (SymbStr::iterator, SymbStr::iterator) noexcept ;	
    bool	isclosed (const_iterator, const_iterator) const noexcept ;
    bool	isclosed () const { return isclosed (cbegin(), cend()) ; }
    SymbStr&	joinends () { resize (joinends (begin(), end())) ; return *this ; }

    static inline bool	dots { false } ;	// print w. symb separators?

    friend ostream& operator<< (ostream&, const SymbStr&) ;
    } ;

namespace Symb					// Symbol namespace
    {
    // Symbol code bits:
    //   0,1:   link direction or fermion type & flavor
    //   2:	reflection (gauge) or conjugation (fermion) flag
    //   3,4,5: symbol type (link, E, Elink, EE, EElink, fermion)
    //   3:	zero step flag
    // masks:
    //   0x03:  link axis
    //   0x04:	reflection flag
    //	 0x07:	direction mask
    //   0x08:	zero step flag
    //   0x38:	symbol type mask
    //   0x3c:	symbol type and reflection mask

    static const symb Rflag { 0x04 } ;	// reflection flag
    static const symb addE  { 0x10 } ;	// link -> Elink, etc
    static const symb addEE { 0x20 } ;	// link -> EElink, etc
    static const symb justE { 0x08 } ;	// first E symbol
    static const symb KinG  { 0x18 } ;	// first EE symbol
    static const symb YMend { 0x28 } ;	// end of pure gauge symbols
    static const symb Nsymb { 0x30 } ;	// end of matter symbols
    static const symb KinF  { 0x30 } ;	// first EE_F symbol
    static const symb TandR { 0x3c } ;	// type and reflection mask
    static const symb EntrG { 0x34 } ;	// gauge entropy symbol
    static const symb EntrF { 0x35 } ;	// fermi entropy symbol
    static const symb Null  { 0x36 } ;	// null ligature flag
    static const symb X     { 0x37 } ;	// invalid ligature flag

    inline bool	isstag	  (symb c)	{ return c & 0x01 ; }
    inline bool	isderiv	  (symb c)	{ return c & 0x01 ; }
    inline int	flav	  (symb c)	{ return c & 0x02 ; }
    inline int  axis	  (symb c)	{ return c & 0x03 ; }
    inline bool isrefl	  (symb c)	{ return c & 0x04 ; }
    inline bool isconj	  (symb c)	{ return c & 0x04 ; }
    inline bool	nostep	  (symb c)	{ return c & 0x08 ; }
    inline int  direction (symb c)	{ return c & 0x07 ; }
    inline int	tnR	  (symb c)	{ return c >> 2  ; }
    inline int	type	  (symb c)	{ return c >> 3  ; }

    inline bool	nonlink	  (symb c)	{ return type(c) != 0 ; }
    inline bool	nongauge  (symb c)	{ return type(c) >  4 ; }
    inline bool isxtra	  (symb c)	{ return type(c) >  5 ; }
    inline bool islink	  (symb c)	{ return type(c) == 0 ; }
    inline bool isE	  (symb c)	{ return type(c) == 1 ; }
    inline bool isElink	  (symb c)	{ return type(c) == 2 ; }
    inline bool isEE	  (symb c)	{ return type(c) == 3 ; }
    inline bool isEElink  (symb c)	{ return type(c) == 4 ; }
    inline bool isferm	  (symb c)	{ return type(c) == 5 ; }

    inline bool isL	  (symb c)	{ return tnR(c) ==  1 ; }
    inline bool isLEl	  (symb c)	{ return tnR(c) ==  3 ; }
    inline bool isLE	  (symb c)	{ return tnR(c) ==  5 ; }
    inline bool isf	  (symb c)	{ return tnR(c) == 10 ; }
    inline bool isF	  (symb c)	{ return tnR(c) == 11 ; }
    inline bool isEE_F	  (symb c)	{ return tnR(c) == 12 ; }
    inline bool isentpy	  (symb c)	{ return tnR(c) == 13 ; }

    inline symb stag	  (symb c)	{ return (c ^ 0x01) ; }
    inline symb refl	  (symb c)	{ return (c ^ 0x04) ; }
    inline symb conj	  (symb c)	{ return type(c) != 1 &&
    						 type(c) != 3 ?
						 refl(c) : c ; }
    inline bool	inclE	  (symb c)	{ return type(c)  > 0 &&
    						 type(c)  < 5 ; }
    inline int	step	  (symb c)	{ return (c & 0x08) ? 0 :
						 (c & 0x04) ? -1 : 1 ; }
    inline bool EorElink  (symb c)	{ return isE(c)  || isElink(c)  ; }
    inline bool EEorEElink(symb c)	{ return isEE(c) || isEElink(c)
							 || isEE_F(c) ; }

    int		to_symb   (char) noexcept ;	// char -> symb
    SymbStr	to_symb   (const string&) ;	// string -> SymbStr
    bool	in_thy    (symb) noexcept ;	// Valid symb?

    static const vector<string>	symbname	// symb -> printable map
	{
	"x",    "y",    "z",    "w",	// 0x00 .. 0x03	(l)ink
	"X",    "Y",    "Z",    "W",	// 0x04 .. 0x07 (L)ink
	"A",    "B",    "C",    "D",	// 0x08 .. 0x0b	E
	"XAx",  "YBy",  "ZCz",  "WDw",	// 0x0c .. 0x0f	LEl
	"Ax",   "By",   "Cz",   "Dw",	// 0x10 .. 0x13	El
	"XA",   "YB",   "ZC",   "WD",	// 0x14 .. 0x17	LE
	"AA",   "BB",   "CC",   "DD",	// 0x18 .. 0x1b	EE
	"XAAx", "YBBy", "ZCCz", "WDDw",	// 0x1c .. 0x1f	LEEl
	"AAx",  "BBy",  "CCz",  "DDw",	// 0x20 .. 0x23	EEl
	"XAA",  "YBB",  "ZCC",  "WDD",	// 0x24 .. 0x27	LEE
	"f",     "g",    "h",    "i",	// 0x28 .. 0x2b	f
	"F",     "G",    "H",    "I",	// 0x2c .. 0x2f	F
	"aa",    "bb",  "cc",   "dd",	// 0x30 .. 0x33	EE_f
	"S",     "s"			// 0x34 .. 0x37	S_g, S_f, X
	} ;

    static const symb ligtable[10][10]		// ligature table
	{
	//  l  L  E LEl El LE EE LEEl EEl LEE
	  { 0, 1, 0, 4, 0, 2, 0,  8,  0,  6 },	// l    0
	  { 1, 0, 5, 0, 3, 0, 9,  0,  7,  0 },	// L    1
	  { 4, 0, 6, 0, 8, 0, X,  X,  X,  X },	// E    2
	  { 0, 5, 0, 7, 0, 9, X,  X,  X,  X },	// LEl  3
	  { 0, 2, 0, 8, 0, 6, X,  X,  X,  X },	// El   4
	  { 3, 0, 9, 0, 7, 0, X,  X,  X,  X },	// LE   5
	  { 8, 0, X, X, X, X, X,  X,  X,  X },	// EE   6
	  { 0, 9, X, X, X, X, X,  X,  X,  X },	// LEEl 7
	  { 0, 6, X, X, X, X, X,  X,  X,  X },	// EEl  8
	  { 7, 0, X, X, X, X, X,  X,  X,  X }	// LEE  9
	//
	//	ligtable maps tnR pair -> codes of 0, 1, X
	//	or tnR code of resulting ligature:
	//
	// 0 :=					-> no join
	// 1 := l.L,      L.l			-> no resulting symb
	// 2 := l.LE,    El.L			-> E
	// 3 := L.El,    LE.l			-> LEl
	// 4 := E.l,      l.LEl			-> El
	// 5 := L.E,    LEl.L			-> LE
	// 6 := E.E,      l.LEE, El.LE,   EEl.L	-> EE
	// 7 := L.EEl,  LEl.LEl, LE.El,   LEE.l	-> ELLl
	// 8 := l.LEEl,   E.El,  El.LEl,   EE.l	-> EEl
	// 9 := L.EE,   LEl.LE,  LE.E,   LEEl.L	-> LEE
	// X := 				-> invalid (excess E's)
	} ;

    inline int ligature (symb a, symb b)	// Combine symbs?
	{
	if (isferm(a) || isferm(b)) return 0 ;
	if (axis(a) == axis(b))
	    if (symb x { ligtable[tnR(a)][tnR(b)] })
		return (x == 1) ? Null : ((x << 2) | axis(a)) ;
	return 0 ;
	}
    } ;

using namespace Symb ;

#endif
