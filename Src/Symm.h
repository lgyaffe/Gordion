#ifndef SYMM_H
#define SYMM_H
#include "List.h"
#include <bitset>

class Symm ;
class SymmSum ;

class SymmTerm : public Term<doub,uint>			// Symmetry index with coeff
    {
    public:
    using	Term<doub,uint>::Term ;
    SymmSum	operator+(SymmTerm) const ;		// Add SymmTerm
    SymmSum	operator-(SymmTerm) const ;		// Subtract SymmTerm
    SymmTerm	operator/(doub z)   const		// Rescale SymmTerm
		    { return SymmTerm { item, coeff / z } ; }
    } ;

using SymmVec = vector<SymmTerm> ;
using Symmmap = hash<Symm,uint,Symmhash,Symm_eq> ;

class SymmSum : public SymmVec				// Sum of SymmTerm's
    {
    public:
    using	SymmVec::SymmVec ;
    SymmSum&	operator+(SymmTerm) ;			// Add to list
    SymmSum&	operator-(SymmTerm) ;			// Subtract from list
    SymmSum	operator*(SymmSum) ;			// Compose SymmSum's
    SymmSum	operator*(SymmTerm) ;			// Multiply by SymmTerm
    } ;

class Symm						// Symmetry transformation
    {
    public:
    array<symb,Nsymb>		map ;			// Symb transform map
    std::bitset<Nsymb+1>	sgn ;			// Symb sign flips
    string			name ;			// Symmetry name
    int				indx = -1 ;		// Symmetry number

    Symm() : map (Symb::symblist) {}			// default constructor

    pair<int,Op> operator()(const Op&) const ;		  // Transform Op
    int		 operator()(const Symm&) const noexcept ; // Compose tranformations
    bool	 operator==(const Symm&) const noexcept ; // Compare transformations
    bool	 isCodd() const { return sgn[Cbit] ; }	  // conjugating symmetry?
    bool	 is_id()  const { return &(*this) == &Symm::list[0] ; }

    static constexpr int Cbit { Nsymb } ;	// conjugating symm bit number
    static Index<Symm>	 list ;			// Lattice (x C) symmetry transforms
    static Symmmap	 trans2indx ;		// transform to index map
#ifdef DELPHI
    static inline vector<uint> inverse ;	// indices of inverse symmetries
#endif

    static void		 symminit  () ;			// Initialize symmetries
    static SymmTerm	 known	  (const string&&) ;	// Return named Symm
    static uint		 knownindx (const string&&) ;	// Return named Symm index

    friend ostream& operator<< (ostream&, const Symm&) ;
    } ;

inline Symmmap		Symm::trans2indx ;		// transform to index map
inline Index<Symm>	Symm::list ;			// Defined Symm's

#endif
