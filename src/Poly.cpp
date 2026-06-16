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

PolyTerm PolyTerm::reindex (ObsList& listb, const ObsList& lista) const
    {
    PolyTerm ans { *this } ;
    for (int k(0) ; k < PSIZ ; ++k)
	{
	if (!ans[k]) break ;

	Obs	o	{ lista (ans[k]) } ;
	auto	tmp	{ listb.is_known (std::move(o)) } ;
	if (tmp.coeff)	{ ans[k] = tmp[0] ; ans.coeff *= tmp.coeff ; }
	else gripe (format("Observable {} missing in list {}!",o.print(),listb.name)) ;
	}
    return ans ;
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

void ObsPoly::push_map (PolyMap& map)			// Add map terms to ObsPoly
    {
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

void PolyRec::add (RecHdr hdr, PolyMap& map)		// Add PolyMap to PolyRec
    {
    ObsPoly tmp { ObsList::obs } ;
    tmp.push_map (map) ;			// copy to ObsPoly for sorting
    hdr.len = tmp.size() ;
    push_back (hdr) ;
    for (auto& term : tmp)
	{
	auto ptr  { cast_to<const Poly*> (&term) } ;
	insert (DataRec::end(), ptr, ptr + Poly::ptermsize) ;
	}
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
