#include "Poly.h"
#include "Gen.h"
#include "Global.h"

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
	    else fatal ("Need higher order monomial!") ;
	    }
	return ans ;
	}
    }

ObsPoly::ObsPoly (ObsList& l, const PolyElem& poly)	// Poly -> ObsPoly
    : list(l)
    {
    if (obslist().neq (global.obs)) abort ("Bad ObsPoly call") ;
    for (auto pptr { poly.begin() } ; pptr < poly.end() ;)
	{
	PolyTerm t { PolyElem::nextterm (pptr) } ;
	if (t.coeff) push_back (t) ;
	}
    shrink_to_fit() ;
    sort() ;
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

void ObsPoly::push_map (PolyMap& map)			// PolyMap -> ObsPoly
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
    const auto& oplist	{ gen.oplist() } ;
    for (auto& t : gen)
	{
	const auto&	op  { oplist[t.item] } ;
	Obs		obs { op } ;
	PolyTerm	tmp { obslist().is_known (std::move(obs)) } ;
	if (tmp.coeff)	add (tmp * t.coeff * gen.coeff) ;
	else return false ;
	}
    return true ;
    }

void PolyRec::add (const ObsPoly& obspoly)		// Add ObsPoly to PolyRec
    {
    uint n(0) ;
    push_back (Element {}) ;
    for (const auto& term : obspoly)
	{
	auto indx { term.item[0] } ;
	push_back (term.coeff) ; ++n ;
	for (int i(0) ; indx && i < PSIZ-1 ; ++i)
	    {
	    auto next { term.item[i+1] } ;
	    if (next) { push_back (-indx) ; ++n ; indx = next ; }
	    else break ;
	    }
	push_back (indx) ; ++n ;
	}
    (&back() - n)->hdr.len = n ;
    }

void PolyRec::add (PolyMap& map)			// Add PolyMap to PolyRec
    {
    ObsPoly obspoly { map.obslist() } ;
    obspoly.push_map (map) ;				// copy to ObsPoly for sort
    add (obspoly) ;
    }

#ifdef DELPHI
void PolyRec::mk_indexlist ()				// Build Obs index list
    {
    indices.clear() ;
    for (auto recptr { begin() } ; recptr != end() ; ++recptr)
	{
	const PolyElem*	ptr { &(*recptr) } ;
	const PolyElem*	end { ptr++->end() } ;
	while (ptr++ < end)
	    {
	    numb indx ;
	    while ((indx = (ptr++)->index) < 0) indices.emplace_back (-indx,-indx) ;
	    if (indx > 0)		        indices.emplace_back ( indx, indx) ;
	    if (indices.do_compress())	        indices.compress() ;
	    }
	}
    indices.compress() ;
    }
#endif

ostream& operator<< (ostream& stream, const PolyIndx& t)	// Print PolyIndx
    {
    char c { '(' } ;
    for (auto k : t) { stream << c << k ; c = ',' ; }
    return stream << ')' ;
    }

void PolyTerm::print (ostream& stream, const ObsList& l) const	// Print PolyTerm
    {
    Print::coeffprt (stream, coeff) ;
    if (!item[0] && Print::is_one (std::abs(coeff))) stream << l(item[0]) ;
    for (int k(0) ; k < PSIZ ; ++k)
	if (item[k]) stream << (k ? " <" : "<") << l(item[k]) << ">" ;
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

ostream& operator<< (ostream& stream, const PolyElem& poly)	// Print packed Poly
    {
    string	sep	{ poly.len() > 6 ? "\n\t" : " " } ;
    int		count	(0) ;

    for (auto pptr { poly.begin() } ; pptr < poly.end() ; ++count)
	{
	stream << sep ;
	PolyTerm t { PolyElem::nextterm (pptr) } ;
	t.print (stream, global.obs) ;
	}
    return stream << (count ? "" : " 0") ;
    }

