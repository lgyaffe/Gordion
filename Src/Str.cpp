#include "Str.h"
#include "Gripe.h"

Str::Str (const string& s)		// Construct from printable string
    {
    auto it { std::find (symbname.begin(), symbname.end(), s) } ;
    if (it == symbname.end())
	{
	resize (s.size()) ;
	auto p { begin() } ;
	for (char c : s)
	    {
	    int x { char_to_symb(c) } ;
	    if (x < 0) gripe (format("Bad symbol character {} in {}", (int) c, s)) ;
	    *p++ = x ;
	    }
	symb x { X } ;
	int  y ;

	auto q { begin() - 1 } ; 
	for (p = begin() ; p < end() ; ++p)
	    {
	    if (x != X && (y = ligature(x,*p)))
		{
		if (y != Null) x = *q = y ;
		else	   x = *(--q) ;
		}
	    else x = *++q = *p ;
	    }
	resize (q + 1 - begin()) ;
	}
    else push_back (it-symbname.begin()) ;
#ifdef DELPHI
    validate() ;
#endif
    }

bool Str::isclosed(const_iterator beg, const_iterator end)	// Closed (sub)string?
    const noexcept
    {
    Coord delta { 0,0,0,0 } ;

    for (auto p { beg } ; p < end ; ++p)
	{
	delta.comp[axis(*p)] += step(*p) ;
	}
    return delta.isclosed (theory.box) ;
    }

int Str::join(const_iterator beg, const_iterator end)	// Append symb range
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

int Str::join(symb x)					// Append single symb
    {
    symb y ;
    if (size() && (y = ligature(back(), x)))
	{
	if (y != Null)	{ back() = y ; return 0 ; }
	else		{ pop_back() ; return -1 ; }
	}
    else		{ push_back (x) ; return 1 ; }
    }

void Str::excise(Str::iterator subb, Str::iterator sube)	// Excise substring
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

int Str::joinends (Str::iterator beg, Str::iterator end)	// Join ends
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

string Str::print() const					// Print -> string
    {
    string buf ;
    for (auto& c : *this) buf += symbname[c] ;
    if (buf.empty()) buf += "1" ;
    return buf ;
    }

string Str::print(const_iterator ptr, const_iterator end) const	// Print substring
    {
    string buf ;
    while (ptr < end) buf += symbname[*ptr++] ;
    if (buf.empty()) buf += "1" ;
    return buf ;
    }

ostream& operator<< (ostream& stream, const Str& str)		// Print -> stream
    {
    if (str.size())
	{
	string	delim { Str::dots ? "." : "" } ;
	string	d { "" } ;
	for (symb c : str)
	    {
	    stream << d << symbname[c] ;
	    d = delim ;
	    }
	}
    else stream << "1" ;
    return stream ;
    }

#ifdef DELPHI
#include "Global.h"
#include "Symm.h"
#include "Blab.h"

void Str::validate () const
    {
    if (size() == 0) return ;
    if (Blab::level(Blab::STR) > 1)
	cout << "Str::validate " << *this << "\n" << flush ;

    int	S	(0) ;
    int	fs	(0) ;
    int	Fs	(0) ;
    int	Es	(0) ;
    int	nf	(0) ;
    int	dim	(0) ;
    int derivs	(0) ;
    bool iseuc	{ theory.euclid } ;

    const char* err { nullptr } ;
    for (const auto c : *this)
	{
	if (islink(c) && dim <= axis(c))	dim = axis(c)+1 ;
	if (isferm(c) && nf  <= flav(c))	nf  = flav(c)/2+1 ;
	if (is_S(c))				++S ;
	if (isferm(c) && !isconj(c))		++fs ;
	if (isferm(c) &&  isconj(c))		++Fs ;
	if (isferm(c) && isderiv(c))		++derivs ;
	if (isE(c)  || isElink(c))		++Es ;
	if (isEE(c) || isEElink(c))		Es += 2 ;
	}
    if      (nf  > theory.nf)	 	err = "has excess fermion flavors" ;
    else if (dim > theory.dim)	 	err = "exceeds lattice dimension" ;
    else if (S && size() != 1)	 	err = "is malformed entropy" ;
    else if (fs + Fs && (fs != Fs))	err = "has bad fermion insertions" ;
    else if (fs > 1)			err = "has excessive fermions" ;
    else if (Fs > 1)			err = "has excessive Fermions" ;
    else if ((fs || Fs) && (Es > 1))	err = "has excessive E's" ;
    else if (Fs && !is_F(front()))	err = "is malformed fermion string" ;
    else if (fs && !is_f(back()))	err = "is malformed fermion string" ;
    else if (iseuc && derivs)		err = "has Grassmann derivative" ;

    if (!err && size() > 1)
	{
	for (auto p = cbegin()+1 ; !err && p < cend() ; ++p)
	    {
	    if (ligature (p[-1],p[0])) err = "is malformed" ;
	    }
	}
    if (err) gripe (format("Bad Str: {} {}", print(), err)) ;
    }

bool Str::nongauge () const				// All symbs gauge?
    {
    for (symb c : *this) if (Symb::nongauge (c)) return true ;
    return false ;
    }

bool Str::alllinks () const				// All link symb's?
    {
    for (symb c : *this) if (!islink(c)) return false ;
    return true ;
    }

int Str::trans (const Symm& symm) noexcept		// Symmetry transform Str
    {
    if (front() == EntrG || front() == EntrF)	return 1 ;
    if (symm.is_id())				return 1 ;

    bool neg	{ false } ;
    uint blab	{ Blab::level(Blab::SYMB) } ;
    if (blab > 0) cout << "Symm " << symm.name << " on " << *this << " -> " ;
    if (symm.isCodd()) std::reverse (begin(), end()) ;
    for (auto& c : *this)
	{
	neg ^= symm.sgn[c] ;
	c = symm.map[c] ;
	}
    if (blab > 0) cout << (neg ? "-" : "") << *this << "\n" << flush ;
    ++global.count().symmtrans ;
    return neg ? -1 : 1 ;
    }
#endif

