#ifndef SYMM_H
#define SYMM_H
#include "Op.h"
#include "Obs.h"
#include <bitset>

class Symm ;
class SymmSum ;

class SymmTerm : public Term<uint>			// Symmetry index with coeff
    {
    public:
    using	Term<uint>::Term ;
    SymmSum	operator+(SymmTerm) const ;		// Add SymmTerm
    SymmSum	operator-(SymmTerm) const ;		// Subtract SymmTerm
    SymmTerm	operator/(doub z)   const		// Rescale SymmTerm
		    { return SymmTerm { item, coeff / z } ; }
    } ;

using SymmVec = vector<SymmTerm> ;

class SymmSum : public SymmVec				// Sum of SymmTerm's
    {
    public:
    using	SymmVec::SymmVec ;
    SymmSum&	operator+(SymmTerm) ;			// Add to list
    SymmSum&	operator-(SymmTerm) ;			// Subtract from list
    SymmSum	operator*(SymmSum) ;			// Compose SymmSum's
    SymmSum	operator*(SymmTerm) ;			// Multiply by SymmTerm
    } ;

struct Symmhash						// Symm hash function 
    {
    std::size_t operator()(const Symm&) const noexcept ;
    } ;

using Symmmap = hash<Symm,uint,Symmhash> ;

class Symm						// Symmetry transformation
    {
    public:
    array<symb,Nsymb>		map ;			// Symb transform map
    std::bitset<Nsymb+1>	sgn ;			// Symb sign flips
    string			name ;			// Symmetry name
    int				indx = -1 ;		// Symmetry number

    Symm() : map(idmap) {}				// default constructor

    string		print() const ;
    pair<int,Op>	operator()(const Op&) const ;		// Transform Op
    int			operator()(const Symm&) const noexcept ;// Compose tranformations
    bool		operator==(const Symm&) const noexcept ;// Compare transformations
    bool		isCodd() const { return sgn[Cbit] ; }	// conjugating symmetry?
    bool		is_id()  const { return &(*this) == &Symm::list[0] ; }

    friend ostream& operator<< (ostream&, const Symm&) ;

    static constexpr int		Cbit = Nsymb ;	// conjugating symm bit number
    static const array<symb,Nsymb>	idmap ;		// Identity transformation map
    static Index<Symm>			list ;		// Lattice (x C) symmetry transforms
    static Symmmap			trans2indx ;	// transform to index map

    static SymmTerm	known (const string&&) ;	// Return named Symm
    static void		symminit () ;			// Initialize symmetries
    } ;

#endif
