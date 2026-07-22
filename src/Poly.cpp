#include "Poly.h"
#include "Gen.h"

std::size_t Polyhash::operator()(const Polyindx& t)	// Polyindx hash function
    const noexcept
    {
    std::hash<uint> hasher ;
    size_t answer { 0 } ;
    for (auto x : t)
	{
	answer ^= hasher(x) + 0x9e3779b9 + (answer << 6) + (answer >> 2) ;
	}
    return answer;
    }

PolyTerm PolyTerm::operator* (const PolyTerm& fac) const // Multiply PolyTerms
    {
    if (!item[0])
	{
	PolyTerm ans { fac } ;
	ans.coeff *= coeff ;
	return ans ;
	}
    else
	{
	PolyTerm ans { *this } ;
	ans.coeff *= fac.coeff ;
	if (fac.item[0]) 
	    {
	    int i(0), j(0) ;
	    while (ans.item[i] && ++i < PSIZ) ;	// count non-trivial terms
	    while (fac.item[j] && ++j < PSIZ) ;	// count non-trivial terms
	    int k { i + j } ;
	    if (k <= PSIZ)
		{
		while (--j >= 0) ans[i+j] = fac[j] ;
		if (k > 1) ans.item.mysort() ;
		}
	    else fatal ("Need higher order poly!") ;
	    }
	return ans ;
	}
    }

void PolyTerm::print (ostream& stream, const ObsList& l) const	// Print PolyTerm
    {
    Print::coeffprt (stream, coeff) ;
    if (!item[0] && Print::is_one (std::abs(coeff))) stream << l(item[0]) ;
    for (int k(0) ; k < PSIZ ; ++k)
	if (item[k]) stream << (k ? " " : "") << l(item[k]) ;
    }

void ObsPoly::sort ()					// Sort ObsPoly terms
    {
    std::sort (begin(), end(),
        [](const PolyTerm &a, const PolyTerm &b) { return a.item < b.item ; }) ;
    }

ObsPoly& ObsPoly::negate ()				// Negate ObsPoly
    {
    for (auto& t : *this) t.coeff *= -1 ;
    return *this ;
    }

void ObsPoly::add (const Poly& poly)			// Add Poly
    {
    if (obslist().neq (ObsList::obs)) abort ("Bad ObsPoly::add call") ;
    for (const auto& t : poly)
	{
	if (t.coeff) emplace_back (t.item, t.coeff) ;
	}
    shrink_to_fit() ;
    sort() ;
    }

void ObsPoly::push_map (PolyMap& map)			// Add PolyMap
    {
    if (obslist().neq (map.obslist())) abort ("Bad ObsPoly::push_map call") ;
    reserve (map.size()) ;
    for (const auto& [indx,coeff] : map)
	{
	if (coeff) emplace_back (indx, coeff) ;
	}
    map.clear() ;
    shrink_to_fit() ;
    sort() ;
    }

bool PolyMap::add_gen (const Gen& gen)			// Add Gen to PolyMap
    {
    const auto& oplist { gen.oplist() } ;
    for (auto& t : gen)
	{
	const auto&	op  { oplist[t.item] } ;
	Obs		obs { op, (ObsType) op.type, op.order, -1 } ;
	PolyTerm	tmp { obslist().is_known (std::move(obs)) } ;
	if (tmp.coeff) add (tmp * t.coeff * gen.coeff) ;
	else return false ;
	}
    return true ;
    }

void PolyRec::add (RecHdr hdr, ObsPoly& poly)	// Add ObsPoly to PolyRec
    {
    hdr.len = poly.size() ;
    push_back (hdr) ;
    for (auto& term : poly)
	{
	auto ptr  { cast_to<const Poly*> (&term) } ;
	insert (DataRec::end(), ptr, ptr + Poly::ptermsize) ;
	}
    }

void PolyRec::add (RecHdr hdr, PolyMap& map)	// Add PolyMap to PolyRec
    {
    ObsPoly poly { ObsList::obs } ;
    poly.push_map (map) ;			// copy to ObsPoly for sorting
    add (hdr, poly) ;
    }

ostream& operator<< (ostream& stream, const Polyindx& t)	// Print Polyindx
    {
    char c { '(' } ;
    for (auto k : t) { stream << c << k ; c = ',' ; }
    return stream << ')' ;
    }

ostream& operator<< (ostream& stream, const ObsPoly& poly)	// Print ObsPoly
    {
    string	sep { poly.size() > 3 ? "\n\t" : " " } ;
    int		count(0) ;
    for (const auto& t : poly)
	{
	if (count++) stream << sep ;
	t.print (stream, poly.obslist()) ;
	}
    return stream << (count ? "" : " 0") ;
    }

ostream& operator<< (ostream& stream, const PolyMap& map)	// Print PolyMap
    {
    string	sep { map.size() > 3 ? "\n\t" : " " } ;
    int		count(0) ;
    for (const auto& [key,coeff] : map)
	{
	PolyTerm t { key, coeff } ;
	if (count++) stream << sep ;
	t.print (stream, map.obslist()) ;
	}
    return stream << (count ? "" : " 0") ;
    }

ostream& operator<< (ostream& stream, const Poly& poly)		// Print Poly
    {
    string	sep  { poly.len > 3 ? "\n\t" : " " } ;
    int		count(0) ;
    for (const auto& t : poly)
	{
	if (count++) stream << sep ;
	t.print (stream, ObsList::obs) ;
	}
    return stream << (count ? "" : " 0") ;
    }
