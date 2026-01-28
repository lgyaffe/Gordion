#include "Theory.h"
#include "Commute.h"
#include "Global.h"
#include "Splice.h"
#include "Gripe.h"
#include "Blab.h"

void Commute::commute_poly (const Gen& gen, const ObsPoly& poly, PolyMap& ans)
    {
    //
    // Evaluate [Gen,ObsPoly]
    //
    uint blab { Blab::blablevel[BLAB::COMMUTE] } ;
    for (auto term : poly)				// loop over ObsPoly terms
	{
	if (term.coeff) commute_term (gen, term, poly.obslist(), ans) ;
	}
    if (blab > 1)	cout << "commute_poly: returning (" << ans << ")\n" ;
    else if (blab)	cout << "commute_poly: returning \n" ;
    }

void Commute::commute_poly (const Gen& gen, const PolyMap& poly, PolyMap& ans)
    {
    //
    // Evaluate [Gen,PolyMap]
    //
    uint blab { Blab::blablevel[BLAB::COMMUTE] } ;
    for (auto [indx,coeff] : poly)			// loop over PolyMap terms
	{
	PolyTerm term { indx, coeff } ;
	if (term.coeff) commute_term (gen, term, poly.obslist(), ans) ;
	}
    if (blab > 1)	cout << "commute_poly: returning (" << ans << ")\n" ;
    else if (blab)	cout << "commute_poly: returning \n" ;
    }

void Commute::commute_term (const Gen& gen, const PolyTerm& term, ObsList& list, PolyMap& ans)
    //
    // Evaluate [Gen,PolyTerm]
    //
    {
    uint blab { Blab::blablevel[BLAB::COMMUTE] } ;
    if (blab) cout << "\n--- commute_term: gen " << gen << ", term " << term << "\n" ;

    ObsList&	anslist	{ ans.obslist() } ;
    bool	newlist	{ list.neq(anslist) } ;
    Polyindx	indx	{ term } ;
    int		sgnprod	{ 1 } ;
    int		sgn[PSIZ] ;

    if (newlist && term[0] && term[1])		// need re-indexing?
	{
	for (int i(0) ; i < PSIZ ; ++i)
	    {
	    if (!term[i]) continue ;
	    PolyTerm reindx { anslist.catalog (list(term[i])) } ;
	    sgnprod *= (sgn[i] = reindx.coeff) ;
	    indx[i]  = reindx[0] ;
	    }
	}
    else for (int i(0) ; i < PSIZ ; ++i) sgn[i] = 1 ;

    for (int i(0) ; i < PSIZ ; ++i)		// for each Obs in PolyTerm...
	{
	if (!term[i]) continue ;

	Obs b { list(term[i]) } ;
	if (isEE_F(b.front())) b.front() -= 0x18 ;

	if (b.is_Entropy())			// handle Entropy special case
	    {
	    for (const PolyTerm& t : gen.reduction)
		{
		ObsList& redulist { gen.reduction.obslist() } ;
		const Obs& oa { redulist(t[0]) } ;
		const Obs& ob { redulist(t[1]) } ;
		if (blab > 2) cout << " oa " << oa.print()
				   << " ob " << ob.print() << "\n" ;
		PolyTerm result { anslist.catalog(oa,ob) } ;
		if (result.coeff)
		    ans.add (result * (gen.coeff * term.coeff * t.coeff)) ;
		}
	    if (blab > 1) cout << "commute_term: returning (" << ans << ")\n" ;
	    }
	else				// collect non-differentiated terms
	    {
	    PolyTerm factor { Polyindx(), 1 } ;
	    for (int j(0) ; j < PSIZ-1 ; ++j)
		{
		factor[j] = j < i ? indx[j] : indx[j+1] ;
		}
	    for (auto opterm : gen)			// for each term of Gen
		{
		const Op& a { Op::list[opterm.item] } ;

		factor.coeff  = gen.coeff * opterm.coeff * term.coeff ;
		factor.coeff *= sgn[i] * sgnprod ;

		do_commute (a, b, factor, anslist, ans) ;	// differentiate
		}
	    if (gen.is_Eloop() && b.is_EEloop())	// [Eloop,EEloop] extras
		{
		if (i || factor[0]) fatal ("EEloop * bogus extra stuff") ;

		do_split (gen.coeff * term.coeff/4, gen.reduction, b, anslist, ans) ;
		}
	    }
	}
    }

void Commute::do_commute (const Op& a, const Obs& b, PolyTerm factor, ObsList& list, PolyMap& ans)
    //
    // Evaluate [Op,Obs] * PolyTerm factor
    //
    {
    if (b.is_Entropy()) return ;

    if (a.is_Eloop())				// a's E differentiates
	{
	do_commuteA (a, b, factor, list, ans) ;
	}
    if (b.is_Eloop())				// b's E differentiates
	{
	do_commuteB (a, b, factor, list, ans) ;
	}
    else if (b.is_EEloop())			// b's front E(s) differentiate
	{
	do_commuteC (a, b, factor, list, ans) ;
	}
    if (b.is_EEloop() || b.is_Efermion())	// b's midE differentiates
	{
	do_commuteD (a, b, factor, list, ans) ;
	}
    if (a.is_Fermion() && b.is_fermi())		// fermions anticommute
	{
	do_commuteE (a, b, factor, list, ans) ;
	}
    ++global.count().commutes ;
    }

void Commute::do_commuteA (const Op& a, const Obs& b, PolyTerm factor, ObsList& list, PolyMap& ans)
    //
    // Evaluate [Eloop Op a, Obs b] * PolyTerm factor --- a[0] differentiates
    //
    {
    uint	blab { Blab::blablevel[BLAB::COMMUTE] } ;
    symb	x { a.front() } ;
    int		dir { axis(x) } ;
    int		tnRa { tnR(x) } ;
    int		bsize ( b.size() ) ;
    short	corder { cordsum (a,b) } ;

    if (blab > 1) cout << "do_commuteA " << a << ", " << b << "\n" ;

    SymbStr buf ;
    buf.reserve (a.size() + bsize + 1) ;

    for (int i(0) ; i < bsize ; ++i)
	{
	symb y { b[i] } ;
	if (dir != axis(y) || isferm(y)) continue ;

	int 		tnRb = tnR(y) ;
	const Splice	*sp = splicetbl[tnRa][tnRb] ;

	for (int k(0) ; k < 2 ; ++k)
	    {
	    if (int coef { sp[k].coeff } )
		{
		symb 	x1 { sp[k].prefix } ;
		symb 	x2 { sp[k].suffix } ;

		buf.clear() ;
		if (i == 0 && b.EorEEloop() && x2 > 1)
		    {
		    // suffix E at front
					buf.join((x2 << 2) | dir) ;
		    if (bsize > 1)	buf.join(b.cbegin()+1, b.cend()) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(a.cbegin()+1, a.cend()) ;
		    }
		else
		    {
		    // b[0] or prefix at front
		    if (i)		buf.join(b.cbegin(), b.cbegin()+i) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(a.cbegin()+1, a.cend()) ;
		    if (x2 >= 0)	buf.join((x2 << 2) | dir) ;
		    if (bsize > i+1)	buf.join(b.cbegin()+i+1, b.cend()) ;
		    }
		if (blab > 1) cout << "commute A: [" << a << ", " << b << "] i "
				   << i << ": " << coef << " " << buf << "\n" ;

		Obs o { buf.joinends(), b.type, corder, -1 } ;
		if (blab > 1) cout << "commute A -> " << coef << " " << o << "\n" ;
		PolyTerm term { list.assess(o) } ;
		if (factor.coeff && term.coeff) ans.add (term * factor * coef) ;
		}
	    }
	}
    }

void Commute::do_commuteB (const Op& a, const Obs& b, PolyTerm factor, ObsList& list, PolyMap& ans)
    //
    // Evaluate -[Eloop Obs b, Op a] * PolyTerm factor --- b[0] differentiates
    //
    {
    uint	blab { Blab::blablevel[BLAB::COMMUTE] } ;
    symb	x { b.front() } ;
    int		dir { axis(x) } ;
    int		tnRb { tnR(x) } ;
    int		asize ( a.size() ) ;
    short	corder { cordsum (a,b) } ;

    if (blab > 1) cout << "do_commuteB " << a << ", " << b << "\n" ;

    SymbStr buf ;
    buf.reserve (asize + b.size() + 1) ;
    
    for (int i(0) ; i < asize ; ++i)
	{
	symb y { a[i] } ;
	if (dir != axis(y) || isferm(y)) continue ;

	int 		tnRa = tnR(y) ;
	const Splice	*sp = splicetbl[tnRb][tnRa] ;

	for (int k(0) ; k < 2 ; ++k)
	    {
	    if (int coef { -sp[k].coeff })
		{
		symb 	x1 { sp[k].prefix } ;
		symb 	x2 { sp[k].suffix } ;

		buf.clear() ;
		if (i == 0 && a.is_Eloop() && x2 > 1)
		    {
		    // suffix E at front
					buf.join((x2 << 2) | dir) ;
		    if (asize > 1)	buf.join(a.cbegin()+1, a.cend()) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(b.cbegin()+1, b.cend()) ;
		    }
		else
		    {
		    // a[0] or prefix at front
		    if (i)		buf.join(a.cbegin(), a.cbegin()+i) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(b.cbegin()+1, b.cend()) ;
		    if (x2 >= 0)	buf.join((x2 << 2) | dir) ;
		    if (asize > i+1)	buf.join(a.cbegin()+i+1, a.cend()) ;
		    }
		if (blab > 1) cout << "commute B: -[" << a << ", " << b << "] i "
				   << i << ": " << coef << " " << buf << "\n" ;

		Obs o { buf.joinends(), (ObsType) a.type, corder, -1 } ;
		if (blab > 1) cout << "commute B -> " << coef << " " << o << "\n" ;
		PolyTerm term { list.assess(o) } ;
		if (factor.coeff && term.coeff) ans.add (term * factor * coef) ;
		}
	    }
	}
    }

void Commute::do_commuteC (const Op& a, const Obs& b, PolyTerm factor, ObsList& list, PolyMap& ans)
    //
    // Evaluate -[EEloop Obs b, Op a] * PolyTerm factor --- b[0] differentiates
    //
    {
    uint	blab { Blab::blablevel[BLAB::COMMUTE] } ;
    symb	x { b.front() } ;
    int		dir { axis(x) } ;
    int		tnRb { tnR(x) } ;
    short	midE { b.middleE() } ;
    int		asize ( a.size() ) ;
    int		bsize ( b.size() ) ;
    short	corder { cordsum (a,b) } ;

    if (blab > 1) cout << "do_commuteC " << a << ", " << b << "\n" ;

    SymbStr buf ;
    buf.reserve (asize + bsize + 1) ;
    
    for (int i(0) ; i < asize ; ++i)		// b[0] differentiates
	{
	symb y { a[i] } ;
	if (isferm(y) || dir != axis(y)) continue ;

	int 		tnRa = tnR(y) ;
	const Splice	*sp = splicetbl[tnRb][tnRa] ;

	for (int k(0) ; k < 2 ; ++k)
	    {
	    if (int coef { -sp[k].coeff })
		{
		symb 	x1 { sp[k].prefix } ;
		symb 	x2 { sp[k].suffix } ;

		buf.clear() ;
		if (a.is_Fermion() || (a.is_Eloop() && i > 0))	// a[0] at front
		    {
		    if (i)		buf.join(a.cbegin(), a.cbegin()+i) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(b.cbegin()+1, b.cend()) ;
		    if (x2 >= 0)	buf.join((x2 << 2) | dir) ;
		    if (asize > i+1)	buf.join(a.cbegin()+i+1, a.cend()) ;
		    }
		else if (x2 > 1)				// suffix E at front
		    {
					buf.join((x2 << 2) | dir) ;
		    if (asize > i+1)	buf.join(a.cbegin()+i+1, a.cend()) ;
		    if (i)		buf.join(a.cbegin(), a.cbegin()+i) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(b.cbegin()+1, b.cend()) ;
		    }
		else if (midE && a.is_Loop())		// b[midE] at front
		    {
					buf.join(b.cbegin()+midE, b.cend()) ;
		    if (x2 >= 0)	buf.join((x2 << 2) | dir) ;
		    if (asize > i+1)	buf.join(a.cbegin()+i+1, a.cend()) ;
		    if (i)		buf.join(a.cbegin(), a.cbegin()+i) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
		    if (midE > 1)	buf.join(b.cbegin()+1, b.cbegin()+midE) ;
		    }
		else 						// prefix at front
		    {
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(b.cbegin()+1, b.cend()) ;
		    if (x2 >= 0)	buf.join((x2 << 2) | dir) ;
		    if (asize > i+1)	buf.join(a.cbegin()+i+1, a.cend()) ;
		    if (i)		buf.join(a.cbegin(), a.cbegin()+i) ;
		    }
		if (blab > 1) cout << "commute C: -[" << b << ", " << a << "] i "
				   << i << " k " << k << ": " << coef << " " << buf << "\n" ;

		Obs o { buf.joinends(), commtype (a, b), corder, -1 } ;
		if (blab > 1) cout << "commute C -> " << coef << " " << o << "\n" ;
		PolyTerm term { list.assess(o) } ;
		if (factor.coeff && term.coeff)
		    {
		    if (blab > 1) cout << "commute C => " << " term " << term
				       << " factor " << factor
				       << " coef " << coef << "\n" ;
		    ans.add (term * factor * coef) ;
		    }
		}
	    }
	}
    }

void Commute::do_commuteD (const Op& a, const Obs& b, PolyTerm factor, ObsList& list, PolyMap& ans)
    //
    // Evaluate -[EEloop/Efermion Obs b, Op a] * PolyTerm factor --- b[midE] differentiates
    //
    {
    short midE { b.middleE() } ;
    if (!midE) return ;

    uint	blab { Blab::blablevel[BLAB::COMMUTE] } ;
    symb	x { b.front() } ;
    int		dir { axis(x) } ;
    int		tnRb { tnR(x) } ;
    int		asize ( a.size() ) ;
    int		bsize ( b.size() ) ;
    short	corder { cordsum (a,b) } ;

    if (blab > 1) cout << "do_commuteD " << a << ", " << b << "\n" ;

    SymbStr buf ;
    buf.reserve (asize + bsize + 1) ;
    
    x    = b[midE] ;
    dir  = axis(x) ;
    tnRb = tnR(x) ;
    
    for (int i(0) ; i < asize ; ++i)
	{
	symb y = a[i] ;
	if (isferm(y) || dir != axis(y)) continue ;

	int 		tnRa = tnR(y) ;
	const Splice	*sp = splicetbl[tnRb][tnRa] ;

	for (int k(0) ; k < 2 ; ++k)
	    {
	    doub coef ( -sp[k].coeff ) ;
	    if (coef)
		{
		symb 	x1 { sp[k].prefix } ;
		symb 	x2 { sp[k].suffix } ;
		int	len (0) ;

		buf.clear() ;
		if (a.type != OpType::Fermion)		// b[0] at front
		    {
					len  += buf.join(b.cbegin(), b.cbegin() + midE) ;
		    if (x2 >= 0)	len  += buf.join((x2 << 2) | dir) ;
		    if (i < asize-1)	len  += buf.join(a.cbegin()+i+1, a.cend()) ;
		    if (i)			buf.join(a.cbegin(), a.cbegin()+i) ;
		    if (x1 >= 0)		buf.join((x1 << 2) | dir) ;
		    if (bsize > midE+1)		buf.join(b.cbegin() + midE+1, b.cend()) ;
		    }
		else					// a[0] at front
		    {
		    if (i)		len  += buf.join(a.cbegin(), a.cbegin()+i) ;
		    if (x1 >= 0)   	len  += buf.join((x1 << 2) | dir) ;
		    if (bsize > midE+1)	len  += buf.join(b.cbegin() + midE+1, b.cend()) ;
						buf.join(b.cbegin(), b.cbegin() + midE) ;
		    if (x2 >= 0)		buf.join((x2 << 2) | dir) ;
		    if (i < asize-1)		buf.join(a.cbegin()+i+1, a.cend()) ;
		    }
		if (blab > 1) cout << "commute D: -[" << b << ", " << a << "] i "
				   << i << ": " << coef << " " << buf << "\n" ;

		if (b.is_Efermion() && a.is_Fermion())
		    {
		    bool  astag	( isstag (a.front()) ^ isstag (a.back()) ) ;
		    bool  bstag	( isstag (b.front()) ^ isstag (b.back()) ) ;
		    Obs   oa	{ buf.begin(), buf.begin()+len, ObsType::Fermion, -1, -1 } ;
		    Obs   ob	{ buf.begin()+len, buf.end(),   ObsType::Fermion, -1, -1 } ;

		    if (oa.isclosed() && ob.isclosed() && astag == bstag)
			{
			int sgn2 { 1 } ;
			if (isstag (a.back())  && (a.size() % 2)) sgn2 *= -1 ;
			if (isstag (b.front()) && (a.size() % 2)) sgn2 *= -1 ;

			auto p	{ buf.begin() + 1 } ;
			auto q	{ buf.begin() + len + 1 } ;
			int  l1	{ buf.joinends (p, q-2) } ;
			int  l2	{ buf.joinends (q, buf.end()-1) } ;
			Obs  oc	{ p, p+l1, ObsType::Loop } ;
			Obs  od	{ q, q+l2, ObsType::Loop } ;
			if (blab > 1) cout << "commute D3 -> " << coef/4 * sgn2 << " "
					   << oc << " " << od << "\n" ;
			PolyTerm terms3 { list.assess(oc) * list.assess(od) } ;
			if (factor.coeff && terms3.coeff)
			    {
			    ans.add (terms3 * factor * coef * sgn2 * 0.25) ;
			    }
			}
		    int   sgn1	{ -1 } ;
		    if (!list.canonicalize || oa.is_coord() && ob.is_coord())
			{
			if (blab > 1) cout << "commute D2 -> " << coef << " "
					   << oa << " " << ob << "\n" ;
			PolyTerm terms1 { list.assess(oa) * list.assess(ob) } ;
			if (factor.coeff && terms1.coeff)
			    {
			    ans.add (terms1 * factor * sgn1 * coef) ;
			    }
			}
		    ob.front() = stag (ob.front()) ;
		    oa.front() = stag (oa.front()) ;
		    if (b.Esublat() ^ !(i % 2) ^ isL(y)) sgn1 *= -1 ;
		     
		    if (!list.canonicalize || oa.is_coord() && ob.is_coord())
			{
			if (blab > 1) cout << "commute D2 -> " << coef << " "
					   << oa << " " << ob << "\n" ;
			PolyTerm terms2 { list.assess(oa) * list.assess(ob) } ;
			if (factor.coeff && terms2.coeff)
			    {
			    ans.add (terms2 * factor * sgn1 * coef) ;
			    }
			}
		    }
		else 
		    {
		    Obs o { buf.joinends(), commtype (a, b), corder, -1 } ;
		    if (blab > 1) cout << "commute D1 -> " << coef << " " << o << "\n" ;
		    PolyTerm term { list.assess(o) } ;
		    if (blab > 1) cout << "commute D1 => " << " term " << term
				       << " factor " << factor
				       << " coef " << coef << "\n" ;
		    if (factor.coeff && term.coeff) ans.add (term * factor * coef) ;
		    }
		}
	    }
	}
    }

void Commute::do_commuteE (const Op& a, const Obs& b, PolyTerm factor, ObsList& list, PolyMap& ans)
    //
    // Evaluate [bilinear Op a, bilinear Obs b] * PolyTerm factor --- endpoints anticommute
    //
    {
    uint	blab	{ Blab::blablevel[BLAB::COMMUTE] } ;
    short	corder	{ cordsum (a,b) } ;
    bool	iseuc	{ theory.euclid } ;

    if (blab > 1) cout << "do_commuteE " << a << ", " << b << "\n" ;

    SymbStr buf ;
    buf.reserve (a.size() + b.size() + 1) ;

    symb x { a.back() } ;
    symb y { b.front() } ;
    if (flav(x) == flav(y) && (!iseuc || isderiv(x)))
	{
	int sgn { 1 } ;
	buf.assign (a.cbegin(), a.cend()-1) ;
	buf.join (b.cbegin()+1, b.cend()) ;
	if (!iseuc && isstag(x) ^ isstag(y))
	    {
	    if (a.size() % 2) sgn *= -1 ;
	    buf.front() = stag(buf.front()) ;
	    }
	if (blab > 1) cout << "commute E1: [" << a << ", " << b << "] -> "
			   << sgn << " " << buf << "\n" ;
	Obs o { buf, b.type, corder, -1 } ;
	if (blab > 1) cout << "commute E1 -> " << sgn << " " << o << "\n" ;
	PolyTerm result { list.assess(o) } ;
	if (factor.coeff && result.coeff) ans.add (result * factor * sgn) ;
	}

    x = b.back() ;
    y = a.front() ;
    if (flav(x) == flav(y) && (!iseuc || isderiv(y)))
	{
	int sgn { -1 } ;
	buf.assign (b.cbegin(), b.cend()-1) ;
	buf.join (a.cbegin()+1, a.cend()) ;
	if (!iseuc && isstag(x) ^ isstag(y))
	    {
	    if (b.oddlen()) sgn *= -1 ;
	    buf.front() = stag(buf.front()) ;
	    }
	if (blab > 1) cout << "commute E2: [" << a << ", " << b << "] -> "
				<< sgn << " " << buf << "\n" ;
	Obs o { buf, b.type, corder, -1 } ;
	if (blab > 1) cout << "commute E2 -> " << sgn << " " << o << "\n" ;
	PolyTerm result { list.assess(o) } ;
	if (factor.coeff && result.coeff) ans.add (result * factor * sgn) ;
	}
    }

void Commute::do_split (doub coeff, const ObsPoly& redu, const Obs& b, ObsList& list, PolyMap& ans)
    //
    // Evaluate coeff * [Eloop Gen::reduction, EEloop Obs b] -> cubic product of loops
    //
    {
    uint		blab { Blab::blablevel[BLAB::COMMUTE] } ;
    static PolyTerm	one  { Polyindx(), 1 } ;
    ObsList		tmplist { "SplitTemp" } ;

    if (redu.size())
	{
	if (blab > 1)
	    {
	    cout << "--- do_split: list " << list.name << " [" << redu << ", " << b << "] "
		 << " coeff " << coeff << "\n" ;
	    }
	ObsList& redulist { redu.obslist() } ;

	for (const auto& t : redu)
	    {
	    for (int k(0) ; k < 2 ; ++k)
		{
		if (!t[k]) continue ;

		const Op  a { redulist(t[k]) } ;	// convert Obs -> Op
		PolyTerm  factor { t[1-k], t.coeff * coeff } ;

		if (factor[0] && list.neq (redulist))
		    factor[0] = list.catalog(redulist(factor[0]))[0] ;

		PolyMap tmpmap { tmplist } ;
		do_commute (a, b, one, tmplist, tmpmap) ;
		for (const auto& [indx,coeff] : tmpmap)
		    {
		    if (coeff)
			{
			const Obs& obs { tmplist(indx[0]) } ;
			if (blab > 1)
			    {
			    cout << "reducing " << obs << " * "
				 << factor * coeff << "\n" ;
			    }
			do_inner (obs, factor * coeff, list, ans) ;
			}
		    }
		}
	    }
	}
    if (blab > 1) cout << "do_split returning\n" ;
    }

void Commute::do_inner (const Obs& a, const PolyTerm factor, ObsList& list, PolyMap& ans)
    //
    // Inner commutator: (Eloop -> tr([E,loop])) * PolyTerm factor
    //
    {
    uint blab { Blab::blablevel[BLAB::COMMUTE] } ;
    if (blab > 1) cout << "--- do_inner: list " << list.name << " "
		       << a << " factor " << factor << "\n" ;

    SymbStr buf ;
    buf.reserve (a.size() + 1) ;
    buf.assign  (a.cbegin(), a.cend()) ;

    symb	x   { buf.front() } ;
    auto	p1  { buf.begin() } ;
    Coord 	delta { 0,0,0,0 } ;
    int		dir1  { axis(x) } ;
    bool	skip  { false } ;

    switch (tnR(x))
	{
	case 2: ++p1 ; break ;						// E
	case 3: *p1-=0x0c ; buf.push_back (x-0x08) ; skip=true ; break ;// LEl
	case 4: *p1-=addE ; break ;					// El
	case 5: ++p1 ; buf.push_back (x-addE) ; break ;			// LE
	default: fatal ("do_inner: inconsistency") ;
	}

    auto p3 { buf.end() } ;
    for (auto p2 {p1} ; p2 < p3 ; ++p2)
	{
	symb	y { *p2 } ;
	int	sgn  { step(y) } ;
	int	dir2 { axis(y) } ;

	if (dir1 != dir2) { delta.comp[dir2] += sgn ; continue ; }

	auto b2 {p2} ;
	if (sgn < 0) { --delta.comp[dir2] ; ++b2 ; }

	if (blab == 5)
	    cout << "buf " << buf.print()
		 << " p1 " << p1 - buf.begin()
		 << " p2 " << p2 - buf.begin()
		 << " y " << symbname[y]
		 << " sgn " << sgn
		 << " p3 " << p3 - buf.begin()
		 << " skip " << skip
		 << " dir1==dir2 " << (dir1 == dir2)
		 << " delta " << delta.comp[0] << "," << delta.comp[1] << "\n" ;

	if (delta.isclosed (theory.box))
	    {
	    auto b1 {p1} ;
	    auto e1 {b2-1} ;
	    auto e2 {p3-1} ;
	    if (!skip || (e1 >= b1 && e2 >= b2))
		{
		while (e1 > b1 && ligature(*e1,*b1) == Null) { ++b1 ; --e1 ; }
		while (e2 > b2 && ligature(*e2,*b2) == Null) { ++b2 ; --e2 ; }

		short corda (-1), cordb (-1) ;
		if (b1 == ++e1) { corda = 0 ; cordb = a.corder ; }
		if (b2 == ++e2) { cordb = 0 ; corda = a.corder ; }
		Obs oa { b1, e1, ObsType::Loop, corda, -1 } ;
		Obs ob { b2, e2, ObsType::Loop, cordb, -1 } ;
		if (blab > 2) cout << "do_inner -> " << sgn << " "
				   << oa << " * " << ob << "\n" ;
		PolyTerm result { list.assess(oa) * list.assess(ob) } ;
		if (blab > 2) cout << "do_inner => " << result << " * "
				   << factor * sgn << "\n" << flush ;
		if (result.coeff && factor.coeff) ans.add (result * factor * sgn) ;
		}
	    }
	if (sgn > 0) ++delta.comp[dir2] ;
	}
    if (blab > 2) cout << ".\n" << flush ;
    }

Gen& Commute::commute_gen (const Gen& gen1, const Gen& gen2, Gen& ans)
    //
    // Evaluate [Gen, Gen] -> Gen
    //
    {
    uint blab { Blab::blablevel[BLAB::COMMUTE] } ;
    if (blab) cout << "\n--- commute_gen: " << gen1 << "," << gen2
		   << ", ans " << ans << "\n" ;

    if (!ans.size())
	{
	ans.order = gen1.order + gen2.order ;
	ans.coeff = gen1.coeff * gen2.coeff ;
	ans.T_odd = gen1.T_odd ^ gen2.T_odd ;
	ans.imag  = gen1.imag  ^ gen2.imag  ;
	}
    doub scale { gen1.coeff * gen2.coeff / ans.coeff } ;
    for (auto t1 : gen1)
	{
	const Op a { Op::list[t1.item] } ;	// N.B. non-ref
	for (auto t2 : gen2)
	    {
	    const Op b { Op::list[t2.item] } ;	// N.B. non-ref
	    op_commute (scale * t1.coeff * t2.coeff, a, b, ans) ;
	    }
	}
    ans.collect() ;
    ans.inner_commute() ;
    if (blab > 1)  cout << "commute_gen: returning " << ans << "\n" ;
    else if (blab) cout << "commute_gen: returning \n" ;
    return ans ;
    }

void Commute::op_commute (doub coeff, const Op& a, const Op& b, Gen& ans)
    //
    // Evaluate coeff * [Op, Op] -> Gen
    //
    {
    if (!coeff) return ;

    uint blab { Blab::blablevel[BLAB::COMMUTE] } ;
    if (blab) cout << "op_commute[Op,Op]: a " << a << " b " << b
		   << " coeff " << coeff << "\n" ;

    SymbStr buf ;
    if (a.size() + b.size() >= buf.capacity())
	{
	buf.reserve (a.size() + b.size() + 1) ;
	}
    ++global.count().commutes ;

    if (a.is_Eloop())
	{
	symb	x { a.front() } ;
	int	dir { axis(x) } ;
	int	tnRa { tnR(x) } ;
	
	for (int i(0) ; i < b.size() ; ++i)	// walk around Op b, a[0] differentiates
	    {
	    symb y { b[i] } ;
	    if (dir != axis(y) || isferm(y)) continue ;

	    int 		tnRb = tnR(y) ;
	    const Splice	*sp = splicetbl[tnRa][tnRb] ;

	    for (int k(0) ; k < 2 ; ++k)
		{
		if (!sp[k].coeff) break ;
		buf.clear() ;

		symb 	x1 { sp[k].prefix } ;
		symb 	x2 { sp[k].suffix } ;

		if (i == 0 && b.is_Eloop() && x2 > 1)	// suffix E at front
		    {
					buf.join((x2 << 2) | dir) ;
		    if (b.size() > 1)	buf.join(b.cbegin()+1, b.cend()) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(a.cbegin()+1, a.cend()) ;
		    }
		else						// b[0] or prefix at front
		    {
		    if (i)		buf.join(b.cbegin(), b.cbegin()+i) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(a.cbegin()+1, a.cend()) ;
		    if (x2 >= 0)	buf.join((x2 << 2) | dir) ;
		    if (b.size() > i+1)	buf.join(b.cbegin()+i+1, b.cend()) ;
		    }
		if (blab > 1)
		    {
		    cout << "op_commute A: [" << a << ", " << b << "] -> "
			 << coeff << " * " << (int)sp[k].coeff << " " << buf << "\n" ;
		    }
		Op	op { buf.joinends(), b.type, (short)(a.order + b.order) } ;
		uint	indx { Op::store (op) } ;
		ans.settype (op) ;
		ans.emplace_back ( OpTerm ( indx, coeff * sp[k].coeff ) ) ;
		}
	    }
	}
    if (b.is_Eloop())
	{
	symb	x	{ b.front() } ;
	int	dir	{ axis(x) } ;
	int	tnRb	{ tnR(x) } ;
	
	for (int i(0) ; i < a.size() ; ++i)	// walk around Op a, b[0] differentiates
	    {
	    symb y { a[i] } ;
	    if (dir != axis(y) || isferm(y)) continue ;

	    int 		tnRa = tnR(y) ;
	    const Splice	*sp = splicetbl[tnRb][tnRa] ;

	    for (int k(0) ; k < 2 ; ++k)
		{
		if (!sp[k].coeff) break ;
		buf.clear() ;

		symb 	x1 { sp[k].prefix } ;
		symb 	x2 { sp[k].suffix } ;

		if (i == 0 && a.is_Eloop() && x2 > 1)	// suffix E at front
		    {
					buf.join((x2 << 2) | dir) ;
		    if (a.size() > 1)	buf.join(a.cbegin()+1, a.cend()) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(b.cbegin()+1, b.cend()) ;
		    }
		else						// a[0] or prefix at front
		    {
		    if (i)		buf.join(a.cbegin(), a.cbegin()+i) ;
		    if (x1 >= 0)	buf.join((x1 << 2) | dir) ;
					buf.join(b.cbegin()+1, b.cend()) ;
		    if (x2 >= 0)	buf.join((x2 << 2) | dir) ;
		    if (a.size() > i+1)	buf.join(a.cbegin()+i+1, a.cend()) ;
		    }
		if (blab > 1)
		    {
		    cout << "op_commute B: [" << a << ", " << b << "] -> "
			 << -coeff << " * " << (int)sp[k].coeff << " " << buf << "\n" ;
		    }
		Op	op { buf.joinends(), a.type, (short)(a.order + b.order) } ;
		uint	indx { Op::store (op) } ;
		ans.settype (op) ;
		ans.emplace_back ( OpTerm ( indx, -coeff * sp[k].coeff ) ) ;
		}
	    }
	}
    if (a.is_Fermion() && b.is_Fermion())
	{
	bool iseuc { theory.euclid } ;
	symb x	   { a.back() } ;
	symb y	   { b.front() } ;
	if (flav(x) == flav(y) && (!iseuc || isderiv(x) ^ isderiv(y)))
	    {
	    int	sgn { 1 } ;
	    buf.assign (a.cbegin(), a.cend()-1) ;
	    buf.join (b.cbegin()+1, b.cend()) ;
	    if (!iseuc && isstag(x) ^ isstag(y))
		{
		if (a.size() % 2) sgn *= -1 ;
		buf.front() = stag(buf.front()) ;
		}
	    if (blab > 1)
		{
		cout << "op_commute C1: [" << a << ", " << b << "] -> "
		     << sgn * coeff << " " << buf << "\n" ;
		}
	    Op		op { buf, a.type, (short)(a.order + b.order) } ;
	    uint	indx { Op::store (op) } ;
	    ans.settype (op) ;
	    ans.emplace_back ( OpTerm ( indx, sgn * coeff ) ) ;
	    }

	x = b.back() ;
	y = a.front() ;
	if (flav(x) == flav(y) && (!iseuc || isderiv(x) ^ isderiv(y)))
	    {
	    int	sgn { -1 } ;
	    buf.assign (b.cbegin(), b.cend()-1) ;
	    buf.join (a.cbegin()+1, a.cend()) ;
	    if (!iseuc && isstag(x) ^ isstag(y))
		{
		if (b.size() % 2) sgn *= -1 ;
		buf.front() = stag(buf.front()) ;
		}
	    if (blab > 1)
		{
		cout << "op_commute C2: [" << a << ", " << b << "] -> "
		     << sgn * coeff << " " << buf << "\n" ;
		}
	    Op		op { buf, b.type, (short)(a.order + b.order) } ;
	    uint	indx { Op::store (op) } ;
	    ans.settype (op) ;
	    ans.emplace_back ( OpTerm ( indx, sgn * coeff ) ) ;
	    }
	}
    if (blab > 1) cout << "op_commute[Op,Op] returning\n" ;
    }
