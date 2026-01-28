#include "Assess.h"
#include "Global.h"
#include "Commute.h"
#include "Rep.h"
#include "Blab.h"
#include "Gripe.h"
#if (0)
#include <algorithm>
#include <cfloat>
#endif

thread_local int depth (0) ;		// Used in debugging output

PolyTerm ObsList::assess (Obs& a)	// classify/store/approx/discard Obs?
    {
    ++global.count().assessed ;

    uint blab { Blab::blablevel[BLAB::ASSESS] } ;
    if (blab > 1) cout << "\nassess depth " << ++depth << ": " << name
		   << " frozen " << frozen() << " approx " << approx
		   << " classify " << classify << " " << a << "\n" << flush ;

    doub	sgn  ( canonicalize ? a.canon() : a.findstart() ) ;
    uint	indx { find (a) } ;

    if (blab > 1) cout << "assess depth " << depth << ": "
		       << (sgn < 0 ? "-" : "") << a << "\n" ;

    if (indx != UINT_MAX)			// found in list
	{
	if (indx < UINT_MAX-1)			// do consistency check
	    {
	    const Obs& b { (*this)(indx) } ;
	    if (classify && a.corder >= 0 && a.corder < b.corder) 
		{
		//fatal (format ("Bad corder in (): indx {} {} cord {} -> {}",
		//	name, indx, b.print(), b.corder, a.corder)) ;
		cout << format ("Warning: Bad corder in ():  indx {} {} cord {} -> {}\n",
			name, indx, b.print(), b.corder, a.corder) ;
		}
	    }
	if (blab > 1)
	    cout << "assess " << depth-- << ":: found at #" << indx << "\n" ;
	++global.count().found ;
	return PolyTerm {Polyindx(indx), sgn} ;
	}
    else if (frozen())			// script generation phase, no additions
	{
	if (approx)
	    {
	    PolyTerm approx = a.approximate (*this) ;
	    if (approx.coeff)			// return approx
		{
		if (blab > 1) cout << "assess " << depth-- << ":: returning approx "
				   << approx.coeff << " * " << approx.item << "\n" << flush ;
		return approx * sgn ;
		}
	    else if (blab > 1) cout << "assess " << depth-- << ":: cannot approx\n" ;
	    }
	else if (blab > 1) cout << "assess " << depth-- << ":: absent " << a << "\n" ;
	return PolyTerm (Polyindx(), 0) ;
	}
    else if (classify)			// Obs generation phase, determine xorder
	{
	const auto& maxord { global.maxord() } ;
	if (a.classify (*this))
	    {
	    if (!maxord || a.order() == maxord)
		{
		a.shrink_to_fit() ;
		if (global.maxthread == 1) store(a) ;
		else			   retain(a) ;
		if (blab > 1) cout << "assess " << depth-- << ":: retained/stored "
				   << a << " -> " << name << "\n" ;
		++global.count().stored ;
		return PolyTerm (Polyindx(UINT_MAX-1), sgn) ;
		}
	    else if (a.order() < maxord)
		cout << (format ("Obs {} order {} < maxord {}",a.print(),a.order(),maxord)) << "\n" ;
	    }
	else if (blab > 1) cout << "assess " << depth-- << ":: classify discard\n" ;
	else if (blab) cout << "Classify discard: " << a << "\n" ;
	++global.count().discarded ;
	return PolyTerm (Polyindx(), 0) ;
	}
    else					// Temp list: just store
	{
	a.shrink_to_fit() ;
	uint indx { store(a) } ;
	if (blab > 1) cout << "assess " << depth-- << ":: immediate store " << a
			   << " -> " << name << " #" << indx << "\n" ;
	return PolyTerm (Polyindx(indx), sgn) ;
	}
    }

bool Obs::classify (ObsList& list)	// Determine Obs expectation order
    {
    uint	blab	{ Blab::blablevel[BLAB::ASSESS] } ;
    bool	success	{ false } ;
    short	xmin	(0) ;

    if (has_Es())	xmin = noEbound (list) ;
    if (!bilinear())	xmin = u1bound (xmin) ;
    if (xmin) ++global.count().bounded ;
    if (blab > 1) cout << "classify " << depth << ": " << *this
		       << " xmin " << xmin << "\n" << flush ;

    if (corder < 0 || corder + xmin <= global.maxord())
	//
	// Allow corder < 0 to enable "call classify" to work
	{
	success = has_Es() ? reduce (xmin,list) : factorize (xmin,list) ;
	}
    if (blab > 1) cout << "classify " << depth << ":: " << *this << "\n" << flush ;
    if (success) ++global.count().classified ;
    return success ;
    }

bool Obs::factorize (short xmin, const ObsList& list)	// Factorize no-E Obs
    const
    {
    uint blab { Blab::blablevel[BLAB::ASSESS] } ;
    if (blab > 1) cout << "factorize " << depth << ": " << *this << "\n" << flush ;

    vector<Site>&	sites { sitelist() } ;
    bool		prev  { list.freezeif() } ;
    bool		isFf  { is_Fermion() } ;
    bool		iseuc { theory.euclid } ;
    auto		s { sites.begin() } ;
    auto		b { cbegin() } ;
    auto		e { cend() } ;
    int			intersects (0) ;

    thread_local SymbStr buf ;
    buf.clear() ;
    buf.reserve (size() + 4 * isFf) ;

    if (corder >= 0 && xorder < 0) xorder = corder ;

    while ((s = std::adjacent_find (s, sites.end())) != sites.end())
	{
	auto t { s + 1 } ;				// For each split...
	do  {
	    ++intersects ;
	    PolyTerm fa, fb ;
	    bool pre  { s->pos > 1 } ;
	    bool post { t->pos < size() - 1 } ;

	    if (blab > 2)
		{
		cout << "factorize: trying " << s->pos << "," << t->pos << " "
		     << print(b,b + s->pos) << '|'
		     << print(b + s->pos, b + t->pos) << '|'
		     << print(b + t->pos, e) << "\n" ;
		}
	    if (isFf && (pre || post))
		{
		symb f ( YMend | flav(buf.front()) ) ;
		symb F ( refl (f) ) ;

		if (pre)
		    {
		    buf.assign (b, b + s->pos + 1) ;
		    bool noflip ( iseuc || isstag(buf.front()) ^ (buf.size() % 2) ) ;
		    buf.back() = noflip ? f : stag(f) ;
		    fa = list.is_known (Obs(buf, ObsType::Fermion, -1, -1)) ;
		    if (blab > 3) cout << "is_known " << buf.print() << " -> " << fa << "\n" ;
		    }
		else fa.coeff = 1 ;

		if (post)
		    {
		    buf.assign (b + t->pos - 1, e) ;
		    bool noflip ( iseuc || isstag(buf.back()) ^ (buf.size() % 2) ) ;
		    buf.front() = noflip ? F : stag(F) ;
		    fb = list.is_known (Obs(buf, ObsType::Fermion, -1, -1)) ;
		    if (blab > 3) cout << "is_known " << buf.print() << " -> " << fb << "\n" ;
		    }
		else fb.coeff = 1 ;

		if (fa.coeff && fb.coeff)
		    {
		    buf.assign (b + s->pos-1, b + t->pos + 1) ;
		    buf.front() = iseuc || (buf.size() % 2) ? F : stag(F) ;
		    buf.back()  = f ;
		    PolyTerm fc { list.is_known (Obs(buf, ObsType::Fermion, -1, -1)) } ;
		    if (blab > 3) cout << "is_known " << buf.print() << " -> " << fc << "\n" ;

		    if (fc.coeff)
			{
			short xord0 ( list(fa[0]).xorder ) ;
			short xord1 ( list(fb[0]).xorder ) ;
			short xord2 ( list(fc[0]).xorder ) ;
			if (xorder < 0 || xord0 + xord1 + xord2 < xorder)
			    {
			    if (blab > 1)
				{
				cout << "factorize " << depth << ": "
				     << *this << " -> "
				     << list(fa[0]) << " "
				     << list(fb[0]) << " "
				     << list(fc[0]) << "\n" << flush ;
				}
			    xorder = xord0 + xord1 + xord2 ;
			    if (xorder == xmin) break ;
			    }
			}
		    else if (blab > 2) cout << "factorizeC " << depth << " unknown factor c\n" ;
		    }
		else if (blab > 2) cout << "factorizeAB " << depth << " unknown factors ab\n" ;
		}
	    buf.assign (b + s->pos, b + t->pos) ;
	    buf.joinends () ;
	    fa = list.is_known (Obs(buf, ObsType::Loop, -1, -1)) ;
	    if (blab > 3) cout << "is_known " << buf.print() << " -> " << fa << "\n" ;
	    if (fa.coeff)
		{
		buf.assign (b, b + s->pos) ;
		buf.join   (b + t->pos, e) ;
		buf.joinends() ;
		fb = list.is_known (Obs(buf, type, -1, -1)) ;
		if (blab > 3) cout << "is_known " << buf.print() << " -> " << fb << "\n" ;

		if (fb.coeff)
		    {
		    short xord0 ( list(fa[0]).xorder ) ;
		    short xord1 ( list(fb[0]).xorder ) ;
		    if (xorder < 0 || xord0 + xord1 < xorder)
			{
			if (blab > 1)
			    {
			    cout << "factorize " << depth << ": " << *this << " -> "
				 << list(fa[0]) << " "
				 << list(fb[0]) << "\n" << flush ;
			    }
			xorder = xord0 + xord1 ;
			if (xorder == xmin) break ;
			}
		    else if (blab > 1)
			{
			cout << "factorize " << depth << ": xorder " << xorder
			     << " xord0 " << xord0 << " xord1 " << xord1 << "\n" ;
			}
		    }
		else if (blab > 2) cout << "factorize " << depth << " unknown factors3\n" ;
		}
	    } while (++t != sites.end() && *t == *s) ;

	if (xorder == xmin) break ;
	s = t ;
	}
    list.refreezeif (prev) ;

    if (!intersects)		++global.count().nointersect ;
    else if (known_xord())	++global.count().factored ;

    if (!known_xord() && blab)
	{
	cout << "Can't factor " << *this << "(" << corder << "," << xorder << ")\n" ;
	}
    if (blab > 1) cout << "factorize " << depth << ":: " << *this << "\n" << flush ;
    return known_xord() ;
    }

bool Obs::reduce (short xmin, ObsList& list)	// Commute w. primaries to remove E's
    {
    PolyTerm	one   { Polyindx(), 1 } ;
    bool	prev  { list.freezeif() } ;
    bool	isham { !theory.euclid } ;
    uint	blab  { Blab::blablevel[BLAB::ASSESS] } ;
    if (blab > 1) cout << "reduce " << depth << ": " << *this << "\n" << flush ;

    if (is_EEloop())
	{
	for (const auto& gen : global.info(0).gens.front())
	    {
	    if (gen.order != 2 || !gen.is_Eloop()) continue ;
	    if (blab > 1)
		{
		cout << "reduce " << depth << ": splitting ["
		     << gen.reduction << "," << *this << "]\n" << flush ;
		}
	    PolyMap ans { list } ;
	    Commute::do_split (1.0, gen.reduction, *this, list, ans) ;
	    for (const auto& [indx,coeff] : ans)
		{
		//if (coeff)		// N.B.
		    {
		    const Obs& o1 { list(indx[0]) } ;
		    const Obs& o2 { list(indx[1]) } ;
		    short xord ( o1.xorder + o2.xorder + 2 ) ;
		    if (xorder < 0 || xorder > xord)
			{
			if (blab > 1) cout << "reduceA: " << *this << " xorder " << xorder
			     << " -> " << xord << "\n" ;
			xorder = xord ;
			}
		    if (xorder == xmin) break ;
		    }
		}
	    if (xorder == xmin) break ;
	    }
	}
    if (xorder != xmin)
	{
	for (const auto& op : Op::list)
	    {
	    if (op.is_Loop())
		{
		if (op.order != 2 || !has_Es()) continue ;
		}
	    else if (op.is_Fermion())
		{
		if (op.order != 1 || op.is_coord() || !bilinear()) continue ;
		}
	    else continue ; // op is Eloop()

	    PolyMap ans { list } ;
	    if (isham && op.is_Loop() && bilinear()) front() = stag(front()) ;
	    if (blab > 1)
		{
		cout << "reduce " << depth << ": doing ["
		     << op << "," << *this << "]\n" << flush ;
		}
	    if (op.is_Fermion() && is_Efermion())
		 Commute::do_commuteD (op, *this, one, list, ans) ;
	    else Commute::do_commute  (op, *this, one, list, ans) ;
	    if (isham && op.is_Loop() && bilinear()) front() = stag(front()) ;

	    for (const auto& [indx,coeff] : ans)
		{
		//if (coeff)		// N.B.
		    {
		    short xord { op.order } ;
		    for (auto& k : indx) if (k) xord += list(k).xorder ;
		    if (xorder < 0 || xorder > xord)
			{
			if (blab > 1) cout << "reduceB: " << *this << " xorder " << xorder
					   << " -> " << xord << "\n" ;
			xorder = xord ;
			}
		    if (xorder == xmin) break ;
		    }
		}
	    if (xorder == xmin) break ;
	    }
	}
    list.refreezeif (prev) ;

    if (blab > 1) cout << "reduce " << depth << ":: " << *this
		       << ", xorder " << xorder << "\n" << flush ;

    if (known_xord()) ++global.count().reduced ;
    return known_xord() ;
    }

PolyTerm Obs::approximate (const ObsList& list) const		// Approximate Obs?
    {
    uint blab { Blab::blablevel[BLAB::ASSESS] } ;
    if (blab > 1) cout << "approximate " << depth << ": " << *this << "\n" ;

    thread_local vector<Intersect> intersects ;
    intersects.clear() ;

    vector<Site>&	sites { sitelist() } ;
    auto		s { sites.begin() } ;

    while ((s = std::adjacent_find (s, sites.end())) != sites.end())
	{
	auto t { s + 1 } ;
	do  {					// Score possible split
	    float score ;
	    short deltaEs  ( t->nEs - s->nEs ) ;

	    if (deltaEs % 2) continue ;

	    if (!bilinear())	// Strongly prefer more even factorizations
		{
		float asym ( 2.0 * (t->pos - s->pos) / size() - 1.0 ) ;
		score = 10 * abs(asym) ;
		}
	    else 		// Strongly prefer large subloop
		{
		score = 10 * abs(t->pos - s->pos) / size() ;
		}
	    // Prefer not to cut at E's
	    if (inclE (s->in) || inclE (s->out)) score += 2 ;
	    if (inclE (t->in) || inclE (t->out)) score += 2 ;

	    if (!isferm (s->in) && !isferm (t->out))
		{ // weakly prefer crossings w/o multiple link traversals
		if (direction (s->in)  == direction (t->in))	   score += 1 ;
		if (direction (s->out) == direction (t->out))	   score += 1 ;
		if (direction (s->out) == direction (refl(t->in))) score += 1 ;
		if (direction (t->out) == direction (refl(s->in))) score += 1 ;
		}
	    intersects.push_back (Intersect{score, s->pos, t->pos, deltaEs}) ;
	    //intersects.emplace_back (score, s->pos, t->pos, deltaEs) ;

	    } while (++t != sites.end() && *t == *s) ;
	s = t ;
	}

    if (intersects.size())
	{
	std::sort (intersects.begin(), intersects.end()) ;

	int attempt ( 0 ) ;

	thread_local SymbStr buf ;
	buf.clear() ;
	buf.reserve (size()) ;

	for (auto& split : intersects)
	    {
	    buf.assign ( cbegin(), cend() ) ;
	    auto	p   { buf.begin() + split.pos1 } ;
	    auto	q   { buf.begin() + split.pos2 } ;
	    int		len { buf.joinends (p,q) } ;

	    ObsType	atype { split.nEs ? ObsType::EEloop : ObsType::Loop } ;
	    ObsType	btype { split.nEs ? ObsType::Loop   : type } ;

	    Obs oa { p, p+len, atype } ;
	    buf.excise (p,q) ;
	    Obs ob { buf.begin(), buf.end(), btype } ;
	    PolyTerm factor { list.is_known (std::move(oa),std::move(ob)) } ;

	    if (!factor.coeff)
		{
		if (++attempt < Intersect::maxtry)	continue ;
		else					break ;
		}
	    if (blab > 1) cout << "approx " << depth << " " << *this << " |"
			       << split.pos1 << "," << split.pos2 << "| -> "
			       << factor.coeff << " * " << factor.item << " = "
			       << list(factor[0]) << " "
			       << list(factor[1]) << "\n" ;
	    ++global.count().approxed ;
	    return factor ;
	    }
	if (blab) cout << "can't approximate: " << *this << "\n" ;
	++global.count().noapprox ;
	}
    else
	{
	if (blab) cout << "no intersect: " << *this << "\n" ;
	++global.count().nointersect ;
	}
    if (blab > 1) cout << "approximate " << depth << ":: " << *this << " failed \n" ;
    return PolyTerm (Polyindx(), 0) ;
    }

vector<Site>& Obs::sitelist () const		// Return sorted list of site coords
    {
    auto	p { cbegin() + 1 } ;
    char	c { front() } ;
    Coord	x {0,0,0,0} ;
    short	nEs (0) ;

    thread_local vector<Site> sites ;
    sites.clear() ;
    sites.reserve (size()) ;

    if (!is_fermi())
	{
	sites.push_back (Site{x.point, 0, 0, back(), c}) ;
	//sites.emplace_back (x.point, 0, 0, back(), c) ;
	char dir ( axis(c) ) ;
	x.comp[dir] += step(c) ;
	x.modlen (dir,theory.box) ;
	if (EorElink(c))	nEs++ ;
	else if (EEorEElink(c))	nEs += 2 ;
	}
    for (short i(1) ; p < cend() ; ++p, ++i)
	{
	sites.push_back (Site{x.point, i, nEs, c, *p}) ;
	//sites.emplace_back (x.point, i, nEs, c, *p) ;
	c = *p ;
	char dir ( axis(c) ) ;
	x.comp[dir] += step(c) ;
	x.modlen (dir,theory.box) ;
	if (EorElink(c)) nEs++ ;
	}
    std::stable_sort (sites.begin(), sites.end()) ;
    return sites ;
    }

short Obs::u1bound (short oldmin) const		// Lower bound xorder, gauge Obs only
    {
    uint	blab { Blab::blablevel[BLAB::ASSESS] } ;
    Coord	bbox[2] {{0,0,0,0}, {0,0,0,0}} ;
    Coord	x {0,0,0,0} ;
    short	minord (0) ;

    if (oldmin == SHRT_MAX) return oldmin ;	// skip if already ruled out
    if (blab > 1) cout << "u1bound " << depth << ": " << *this << "\n" ;

    thread_local vector<vector<int>> flux ;

    for (auto c : *this)
	{
	char dir ( axis(c) ) ;
	x.comp[dir] += step(c) ;
	if (x.comp[dir] < bbox[0].comp[dir]) bbox[0].comp[dir] = x.comp[dir] ;
	if (x.comp[dir] > bbox[1].comp[dir]) bbox[1].comp[dir] = x.comp[dir] ;
	}
    for (int i(0) ; i < theory.dim ; ++i)	// Find flux through ij planes
	{
	for (int j(i) ; ++j < theory.dim ;)
	    {
	    int boxi { bbox[1].comp[i] - bbox[0].comp[i] } ;
	    int boxj { bbox[1].comp[j] - bbox[0].comp[j] } ;
	    flux.clear() ;
	    flux.resize (boxi) ;
	    for (auto& f : flux) { f.clear() ; f.resize (boxj) ; }

	    for (int i(0) ; i < theory.dim ; ++i)
		x.comp[i] = -bbox[0].comp[i] ;

	    for (auto c : *this)
		{
		int delta { step(c) } ;
		if (!delta) continue ;

		int dir   { axis(c) } ;
		if (dir == i)
		    {
		    if (delta < 0) x.comp[dir] += delta ;
		    for (int k (x.comp[j]) ; k < boxj ; ++k)
			{
			flux.at(x.comp[i]).at(k) += delta ;
			}
		    if (delta > 0) x.comp[dir] += delta ;
		    }
		else x.comp[dir] += delta ;
		}
	    for (auto& f : flux)		// Sum up [ij] flux
		{
		for (auto v : f) minord += 2 * abs (v) ;
		}
	    }
	}
    if (blab > 1) cout << "u1bound " << depth << ":: " << *this
		       << " -> " << minord << "\n" ;
    return std::max (oldmin, minord) ;
    }

short Obs::noEbound (const ObsList& list) const	// Lower bound xorder, Obs w. E's
    {
    auto	p { cbegin() } ;
    uint	blab { Blab::blablevel[BLAB::ASSESS] } ;
    if (blab > 1) cout << "noEbound " << depth << ": " << *this << "\n" ;

    thread_local SymbStr buf ;
    buf.clear() ;
    buf.reserve (size()) ;

    switch (Symb::type(*p))				// remove E's
	{
	case 1 : ++p ; break ;					// E, LEl
	case 3 : ++p ; break ;					// EE, LEEl
	case 2 : buf.push_back (*p - addE)  ; ++p ; break ;	// El, LE
	case 4 : buf.push_back (*p - addEE) ; ++p ; break ;	// EEl, LEE
	default: buf.push_back (*p)         ; ++p ; break ;	// F
	}
    auto q { p } ;
    while (q < cend() && islink(*q)) ++q ;
    buf.join (p,q) ;
    if (q < cend())
	{
	switch (Symb::type(*q))
	    {
	    case 1 : ++q ; break ;				// E, LEl
	    case 2 : buf.join (*q - addE) ; ++q ; break ;	// El, LE
	    default: break ;
	    }
	buf.join (q, cend()) ;
	}
    buf.joinends() ;
    if (bilinear()) buf.front() = stag (buf.front()) ;

    ObsType	ntype { is_Efermion() ? ObsType::Fermion : ObsType::Loop } ;
    Obs		obs   { buf, ntype, corder, -1 } ; obs.canon() ;
    uint	indx  { list.find (obs) } ;
    short	xord  ( indx == UINT_MAX ? SHRT_MAX : list(indx).xorder ) ;
    if (blab > 1)
	{
	cout << "noEbound " << depth << ":: " << *this << " -> "
	     << buf << " -> " << obs << " indx " << indx
	     << " returning xord " << xord << "\n" ;
	}
    return xord ;
    }
