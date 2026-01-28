#include "Symb.h"
#include "Gripe.h"

bool Symb::in_thy(symb c) noexcept		// Test if valid symb
    {
    if (c >= Nsymb)				return false ;
    if (c >= YMend && flav(c) >= 2*theory.nf)	return false ;
    if (c <  YMend && axis(c) >=  theory.dim)	return false ;
    return true ;
    }

int Symb::to_symb(char c) noexcept		// Map character -> symb
    {
    for (symb x(0) ; x < symbname.size() ; ++x)
	{
	if (symbname[x].size() > 1) continue ;
	if (symbname[x].front() == c) return x ;
	}
    return -1 ;
    }

SymbStr Symb::to_symb(const string& s)		// Map printables -> SymbStr
    {
    SymbStr a ; a.resize (s.size()) ;

    auto p { a.begin() } ;
    for (char c : s)
	{
	int x { to_symb(c) } ;
	if (x < 0) gripe (format("Bad symbol character {} in {}", (int) c, s)) ;
	*p++ = x ;
	}
    symb x { X } ;
    int  y ;

    auto q { a.begin() - 1 } ; 
    for (p = a.begin() ; p < a.end() ; ++p)
	{
	if (x != X && (y = ligature(x,*p)))
	    {
	    if (y != Null) x = *q = y ;
	    else	   x = *(--q) ;
	    }
	else x = *++q = *p ;
	}
    a.resize (q + 1 - a.begin()) ;
    while (a.size() > 1 && (y = ligature(a.back(),a.front())))
	{
	if (y != Null)
	    {
	    a.front() = y ;
	    a.resize(a.size() - 1) ;
	    }
	else
	    {
	    rotate(a.begin(), a.begin() + 1, a.end()) ;
	    a.resize(a.size() - 2) ;
	    }
	}
    return a ;
    }

bool SymbStr::isclosed(const_iterator beg, const_iterator end)	// Closed (sub)string?
    const noexcept
    {
    Coord delta { 0,0,0,0 } ;

    for (auto p { beg } ; p < end ; ++p)
	{
	delta.comp[axis(*p)] += step(*p) ;
	}
    return delta.isclosed (theory.box) ;
    }

int SymbStr::join(const_iterator beg, const_iterator end)	// Append symb range
    {
    int  k(0) ;
    symb y ;
    while (size() && beg < end && (y = ligature(back(), *beg)))
	{
	if (y != Null)	back() = y ;
	else		{ pop_back() ; --k ; }
	++beg ;
	}
    append (beg, end) ;
    return k + end - beg ;
    }

int SymbStr::join(symb x)					// Append single symb
    {
    symb y ;
    if (size() && (y = ligature(back(), x)))
	{
	if (y != Null)	{ back() = y ; return 0 ; }
	else		{ pop_back() ; return -1 ; }
	}
    else		{ push_back (x) ; return 1 ; }
    }

void SymbStr::excise(SymbStr::iterator subb, SymbStr::iterator sube)	// Excise substring
    noexcept
    {
    auto p { subb - 1 } ;
    auto q { sube } ;
    symb y ;
    while (p >= begin() && q < end() && (y = ligature(*p,*q)))
	{
	if (y != Null)	*p = y ;
	else		--p ;
	++q ;
	}
    while (q < end()) *++p = *q++ ;
    resize (joinends (begin(), p+1)) ;
    }

int SymbStr::joinends (SymbStr::iterator beg, SymbStr::iterator end)	// Join ends
    noexcept
    {
    for (symb y ; end - beg > 1 && (y = ligature(*(end-1),*beg)) ;)
	{
	if (y != Null)
	    {
	    *beg = y ; --end ;
	    }
	else
	    {
	    rotate(beg, beg + 1, end) ; end -= 2 ;
	    }
	}
    if (end - beg == 1)
	{
	symb x { *beg } ;
	if ((isE(x) || isEE(x)) && isrefl(x)) *beg = refl(x) ;
	}
    return end - beg ;
    }

string SymbStr::print() const						// Print -> string
    {
    string buf ;
    for (auto& c : *this) buf += symbname[c] ;
    if (buf.empty()) buf += "1" ;
    return buf ;
    }

string SymbStr::print(const_iterator ptr, const_iterator end) const	// Print substring
    {
    string buf ;
    while (ptr < end) buf += symbname[*ptr++] ;
    if (buf.empty()) buf += "1" ;
    return buf ;
    }

ostream& operator<< (ostream& stream, const SymbStr& str)		// Print -> stream
    {
    if (str.size())
	{
	string delim { SymbStr::dots ? "." : "" } ;
	symb c	     { str.front() } ;
	string d     { "" } ;

	for (char c : str)
	    {
	    stream << d << symbname[c] ;
	    d = delim ;
	    }
	}
    else stream << "1" ;
    return stream ;
    }
