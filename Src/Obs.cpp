#include "Obs.h"
#include "Global.h"
#include "Numerics.h"
#include "Op.h"
#include "Symm.h"
#include "Gripe.h"
#include "Blab.h"
#include <regex>

Obs::Obs(const Op& a)					// Convert Op -> Obs
    :
    Str(a), corder(a.order), xorder(-1)
    {
    type = a.is_Fermion() ? ObsType::Fermion :
	   a.is_Eloop()   ? ObsType::Eloop : ObsType::Loop ;
    }

Obs::Obs(const string& s)				// Construct Obs
    :
    Str(s), type(obstype(s)), corder(-1), xorder(-1)
    {
    if (Blab::level(Blab::OBS) > 3)
	{
	cout << "Obs(string&) " << s << "\n" << flush ;
	}
    joinends() ;
    validate() ;
    }

Obs::Obs(const string& s, ObsType t, short cord, short xord)	// Construct Obs
    :
    Str(s), type(t), corder(cord), xorder(xord)
    {
    if (Blab::level(Blab::OBS) > 3)
	{
	cout << "Obs(string&,ObsType,short,short) " << s << "\n" << flush ;
	}
    joinends() ;
    validate() ;
    }

void Obs::validate () const				// Test if valid Obs
    {
    if (size() == 0 && is_Loop()) return ;
    if (Blab::level(Blab::OBS) > 1)
	cout << "Obs::validate " << *this << "\n" << flush ;

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
    else if (fs + Fs > 2)		err = "has excessive fermions" ;
    else if (fs + Es > 2)		err = "has excessive E's" ;
    else if (fs == 0 && !isclosed())	err = "is not closed loop" ;
    else if (fs && !is_F(front()))	err = "is malformed fermion bilinear" ;
    else if (fs && !is_f(back()))	err = "is malformed fermion bilinear" ;
    else if (Es && islink(front()))	err = "is mis-rotated" ;
    else if (iseuc && derivs)		err = "has Grassmann derivative" ;
    else if (fs &&  Es  && type != ObsType::efermion) err = "has wrong type" ;
    else if (fs && !Es  && type != ObsType::Fermion)  err = "has wrong type" ;
    else if (!fs && !Es && type != ObsType::Loop)     err = "has wrong type" ;
    else if (!fs && Es==1 && type != ObsType::Eloop)  err = "has wrong type" ;
    else if (!fs && Es>1  && type != ObsType::EEloop) err = "has wrong type" ;

    if (!err && size() > 1)
	{
	if (ligature (back(),front()))	err = "is malformed" ;
	for (auto p = cbegin()+1 ; !err && p < cend() ; ++p)
	    {
	    if (ligature (p[-1],p[0])) err = "is malformed" ;
	    }
	}
    if (err) gripe (format("Bad Obs: {} {}", print(), err)) ;
    }

short Obs::middleE() const				// Set & return midE location
    {
    if (midE < 0)
	{
	if (!doubleE(front()))
	    {
	    for (auto ptr = cbegin() ; ++ptr < cend() ;)
		{
		if (singleE (*ptr)) return midE = ptr - cbegin() ;
		}
	    fatal (format("Missing midE in obs {}, type {}", print(), (int)type)) ;
	    }
	else midE = 0 ;
	}
    return midE ;
    }

bool Obs::Esublat() const				// E sublattice ?
    {
    short midE { middleE() } ;
    symb  c { (*this)[midE] } ;
    return is_efermion() ? !(midE % 2) ^ (is_LE(c) || is_LEl(c)) : false ;
    }

ObsType Obs::obstype (const string s)			// Determine Obs type
    {
    std::regex			E_s ("[ABCD]") ;
    std::regex			F_s ("[FGHI]") ;
    std::sregex_iterator	end ;
    std::sregex_iterator	Ebeg (s.begin(), s.end(), E_s) ;
    std::sregex_iterator	Fbeg (s.begin(), s.end(), F_s) ;
    int				nE ( std::distance (Ebeg, end) ) ;
    int				nF ( std::distance (Fbeg, end) ) ;

    if      (nE == 0 && nF == 0) return ObsType::Loop ;
    else if (nE == 1 && nF == 0) return ObsType::Eloop ;
    else if (nE == 2 && nF == 0) return ObsType::EEloop ;
    else if (nE == 0 && nF == 1) return ObsType::Fermion ;
    else if (nE == 1 && nF == 1) return ObsType::efermion ;
    else gripe (format("Obs {} is invalid type", s)) ;
    }

int Obs::findstart() noexcept				// Rotate to preferred start
    {
    if (is_Loop() && size() > 0)
	{
	int a(0) ;
	int b(0) ;
	int len ( size() ) ;
	auto s { c_str() } ;

	while (++b < len)
	    {
	    int k(0) ;
	    while (k < len && s[(a+k) % len] == s[(b+k) % len]) ++k ;
	    if (k < len && s[(b+k) % len] < s[(a+k) % len]) a = b ;
	    }
	if (a) rotate (begin(), begin() + a, end()) ;
	}
    else if (is_EEloop() && !doubleE(front()))
	{
	const auto&	blab { Blab::level(Blab::OBS) } ;
	int		a(0), b(middleE()), k(0) ;
	int		len ( size() ) ;
	auto		s { c_str() } ;

	while (k < len && s[(a+k) % len] == s[(b+k) % len]) ++k ;
	if (k < len && s[(b+k) % len] < s[(a+k) % len])
	    {
	    if (blab > 2)
		{
		cout << "Obs::findstart: rotating obs " << *this << " "
		     << b << " places " << flush ;
		}
	    rotate (begin(), begin() + b, end()) ;
	    midE = len - b ;
	    if (blab > 2) cout << " -> " << *this << "\n" ;
	    }
	}
    return 1 ;
    }

int Obs::trans (const Symm& symm, int start) noexcept	// Symmetry transform Obs
    {
    if (front() == EntrG || front() == EntrF)	return 1 ;
    if (!start && symm.is_id())			return 1 ;

    const auto&	blab	{ Blab::level(Blab::SYMM) } ;
    bool	neg	{ false } ;
    if (blab > 1) cout << "Symm " << symm.name << " on " << *this
		       << " at " << start << " -> " ;
    if (symm.isCodd())
	{
	std::reverse (begin(), end()) ;
	start = size() - 1 - start ;
	midE  = -1 ;
	}
    if (start)
	{
	rotate (begin(), begin() + start, end()) ;
	midE = -1 ;
	}

    for (auto& c : *this)
	{
	neg ^= symm.sgn[c] ;
	c = symm.map[c] ;
	}
    if (!theory.euclid && isferm(back()) && isstag(back()))
	{
	neg ^= oddlen() ;
	back()  = stag(back()) ;
	front() = stag(front()) ;
	}
    if (blab > 1) cout << (neg ? "-" : "") << *this << "\n" << flush ;

    ++global.count().symmtrans ;
    return neg ? -1 : 1 ;
    }

ostream& operator<< (ostream& stream, const Obs& obs)	// Print Obs -> stream
    {
    if (Blab::level(Blab::OBS))
	{
	stream  << "(" << obs.corder << "," << obs.xorder ;
	stream << "," << Obs::type_name [(int)obs.type] ;
	stream << ") " ;
	}
    if (obs.imag()) stream << "i " ;
    stream << (Str) obs ;
    return stream ;
    }

ObsStats::ObsStats (const ObsList& list)		 // Construct ObsStats
    {
    maxloop = 0 ;
    doub	maxloopvev { 0.0 } ;
    ulong	lengthsum (0) ;
    auto&	vevs { numerics.vev  } ;
    auto	nvev { vevs.size()   } ;
    auto 	beg  { list.cbegin() } ;
    auto 	end  { list.cend()   } ;

    for (auto ptr { beg } ; ptr < end ; ++ptr)
	{
	int	indx	( ptr - beg ) ;
	auto	obs	{ list(indx) } ;
	int	t	{ (int) obs.type } ;
	short	c	( std::max (obs.corder, (short) 0) ) ;
	short	x	( std::max (obs.xorder, (short) 0) ) ;
	
	if (c >= (*this)[t].size())	(*this)[t].resize (c+1) ;
	if (x >= (*this)[t][c].size())	(*this)[t][c].resize (x+1) ;

	if (c > maxc) maxc = c ;
	if (x > maxx) maxx = x ;

	++(*this)[t][c][x] ;
	lengthsum += obs.size() ;
	if (obs.is_Loop() && indx && indx < nvev)
	    {
	    if (std::abs(vevs[indx]) > maxloopvev)
		{
		maxloopvev = std::abs(vevs[indx]) ;
		maxloop = indx ;
		}
	    }
	}
    avglen = (float) lengthsum / list.size() ;
    }

