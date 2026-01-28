#include "Canon.h"
#include "Symm.h"
#include "Gripe.h"
#include "Blab.h"
#include <cstring>

int Obs::canon()				// Canonicalize observable
    {
    uint blab { Blab::blablevel[BLAB::CANON] } ;
    if (blab > 1) cout << "canon " << *this << "\n" << flush ;
    ++global.count().canons ;

    int len (size()) ;

    if (!len)					// identity loop
	{
	xorder = corder = 0 ;
	if (blab > 1) cout << " -> " << *this << "\n" << flush ;
	return 1 ;
	}
    else if (isxtra(front()))			// entropy, EE_F special cases
	{
	if (isEE_F(front()))			// fermion contribution to EE
	    {
	    short xord { xorder } ;
	    front() -= 0x18 ;			// convert to regular EE
	    canon() ;				// recurse
	    front() += 0x18 ;			// convert back to EE_F
	    xorder = xord ;
	    }
	if (blab > 1) cout << " -> " << *this << "\n" << flush ;
	return 1 ;
	}
    else if (is_Loop())					// ObsType::Loop
	{
	if (len <= loopchunksize)				// lookup in cache
	    {
	    auto search = Canon::cache.find (*this) ;
	    if (search == Canon::cache.end())		// cache miss
		{
		if (blab > 2) cout << "cache miss: " << *this << "\n" << flush ;
		uint minscore { UINT_MAX } ;
		Obs  best { *this } ;

		for (int rot(0) ; rot < len ; ++rot)
		    {
		    for (const auto& symm : Symm::list)
			{
			uint score (0) ;
			Obs  obs { *this } ;
			obs.trans (symm, rot) ;

			for (auto& x : obs) score = (score << loopbits) + x ;
			if (score < minscore)
			    {
			    minscore = score ;
			    best = obs ;
			    }
			}
		    }
		*this = best ;

		if (!Canon::cache.freeze && corder >= 0)
		    {
		    if (!known_xord()) factorize (0, ObsList::obs) ;
		    if (known_xord() && order() <= global.maxord())
			{
			shrink_to_fit() ;
			if (blab) cout << "canon: storing " << *this << "\n" ;
			uint	indx { ObsList::obs.store (*this) } ;
			int	maxrot { theory.dim > 1 ? len : 0 } ;

			for (int rot(0) ; rot < maxrot ; ++rot)	// add cache entries
			    {
			    for (const auto& symm : Symm::list)
				{
				Obs obs { *this } ;
				obs.trans (symm, rot) ;
				auto p { Canon::cache.try_emplace (obs, indx) } ;
				if (p.second && blab > 1)
				    cout << "cached: " << obs << " -> "
					 << *this << "\n" << flush ;
				}
			    }
			}
		    else if (blab > 2) cout << "Canon: cannot cache " << *this << "\n" ;
		    }
		else if (blab > 2) cout << "Canon: cannot cache " << *this << "\n" ;
		if (blab > 1) cout << "canon L1 returning: " << *this << "\n" << flush ;
		++global.count().cachemiss ;
		}
	    else					// cache hit
		{
		*this = ObsList::obs (search->second) ;
		if (blab > 1) cout << "canon L0 returning: " << *this << "\n" << flush ;
		++global.count().cachehits ;
		}
	    }
	else						// compare starting chunks
	    {
	    int		dim { theory.dim } ;
	    auto&	symmset { Canon::symmset } ;
	    uint	minscore { UINT_MAX } ;
	    uint	code(0) ;
	    int		start(0) ;
	    int		end(0) ;
	    auto	p { data() } ;

	    vector<pair<int,int>> trial ;

	    for (; end < loopchunksize ;)		// fill initial loopchunk
		{
		symb c { p[end++] } ;
		int  x { c >= dim ? refl(c) + dim : c } ;
		code = 2 * dim * code + x ;
		}

	    while (true)				// walk around loop
		{
		Node& node { Canon::looptable[code] } ;

		if (node.score < minscore)
		    {
		    trial.resize(0) ;
		    trial.emplace_back (start, node.indx) ;
		    minscore = node.score ;
		    }
		else if (node.score == minscore)
		    {
		    trial.emplace_back (start, node.indx) ;
		    }
		if (++start == len || dim == 1) break ;

		symb c { p[end++ % len] } ;
		int  x { c >= dim ? refl(c) + dim : c } ;
		code = (2 * dim * code + x) % looptblsize() ;
		}

	    if (trial.size() == 1 && symmset[trial[0].second].size() == 1)
		{
		start = trial[0].first ;
		end   = (start + loopchunksize - 1) % len ;
		auto		symmidx { symmset[trial[0].second][0] } ;
		const auto&	symm { Symm::list [symmidx] } ;
		int		beg  { symm.isCodd() ? end : start } ;
		trans (symm, beg) ;
		if (blab > 2)
		    {
		    cout << "loop walk " << symm.name << " at "
			 << beg << " -> " << *this << "\n" << flush ;
		    }
		}
	    else				// compare remaining possibilities
		{
		int bestsymm {-1} ;
		int beststart ;
		Obs best { *this } ;

		for (auto [start,setindx] : trial)
		    {
		    end   = (start + loopchunksize - 1) % len ;

		    for (auto symmidx : symmset[setindx])
			{
			const auto&	symm { Symm::list [symmidx] } ;
			int  		beg  { symm.isCodd() ? end : start } ;
			Obs		obs  { *this } ;
			obs.trans (symm, beg) ;
			if (blab > 2)
			    {
			    cout << "loop try " << symm.name
				 << " at " << beg << " -> " << obs << "\n" ;
			    }
			if (bestsymm < 0 || memcmp (best.data(), obs.data(), size()) > 0)
			    {
			    bestsymm  = symmidx ;
			    beststart = beg ;
			    best = obs ;
			    }
			}
		    }
		if (blab > 2)
		    {
		    cout << "loop walk " << Symm::list[bestsymm].name
			 << " at " << beststart << " -> " << best << "\n" << flush ;
		    }
		*this = best ;
		}
	    if (blab > 1) cout << "canon L2 returning " << *this << "\n" << flush ;
	    }
	return 1 ;
	}
    else						// type != ObsType::Loop
	{
	int	sgn0 { 1 } ;
	int	sgn1 { 1 } ;
	bool	isFf { is_fermi() } ;
	short	mid  ( midEtype() ? middleE() : 0 ) ;

	if (!theory.euclid && isFf && isstag (back()))
	    {
	    front() = stag (front()) ;
	    back()  = stag (back()) ;
	    if (oddlen()) sgn0 *= -1 ;
	    }
	if (theory.dim == 1 && mid)			// use 1D Gauss to move E
	    {
	    auto p { begin() + mid } ;
	    if (blab > 1) cout << "Converting " << *this ;
	    if (is_EEloop())
		{
		mid = midE = 0 ;
		front() += addE ;
		if (isElink(*p)) *p -= addE ;
		else		 excise(p,p+1) ;
		}
	    else if (mid != size() - 2) // is_Efermion
		{
		mid = midE = size() - 2 ;
		(*this)[mid] += addE ;
		if (isElink(*p)) *p -= addE ;
		else		 excise(p,p+1) ;
		}
	    if (blab > 1) cout << " -> " << *this << "\n" << flush ;
	    }

	if (size() - 2 * isFf  <= specchunksize)	// lookup in cache
	    {
	    auto search = Canon::cache.find (*this) ;
	    if (search == Canon::cache.end())		// cache miss
		{
		if (blab > 2) cout << "cache miss: " << *this << "\n" << flush ;

		int	maxtry { (isFf || !mid) ? 1 : 2 } ;
		Obs	best   { *this } ;
		uint	minscore { UINT_MAX } ;

		for (int shot(0) ; shot < maxtry ; ++shot)
		    {
		    for (const auto& symm : Symm::list)
			{
			bool	flip  { isFf && symm.isCodd() } ;
			int	start ( shot ? mid : (flip ? size() - 1 : 0) ) ;
			Obs	obs   { *this } ;
			int	sgn2  { obs.trans (symm, start) } ;
			uint	score (0) ;

			for (auto& x : obs)
			    {
			    if (!isferm(x)) score = (score << specbits) + x ;
			    }
			if (score < minscore)
			    {
			    minscore = score ;
			    best = obs ;
			    sgn1 = sgn2 ;
			    }
			}
		    }
		*this = best ;
		if (mid) mid = best.middleE() ;

		if (!Canon::cache.freeze && corder >= 0)
		    {
		    if (!known_xord())
			{
			if (has_Es())	reduce    (0, ObsList::obs) ;
			else		factorize (0, ObsList::obs) ;
			}
		    if (known_xord() && order() <= global.maxord())
			{
			shrink_to_fit() ;
			if (blab) cout << "canon: storing " << *this << "\n" ;
			uint indx { ObsList::obs.store (*this) } ;

			for (int shot(0) ; shot < maxtry ; ++shot) // add cache entries
			    {
			    for (const auto& symm : Symm::list)
				{
				bool	flip  { isFf && symm.isCodd() } ;
				int	start ( shot ? mid : (flip ? size()-1 : 0) ) ;
				Obs	obs   { *this } ;
				int	sgn2  { obs.trans (symm, start) } ;
				int	code  ( sgn1 * sgn2 < 0 ? -indx : indx ) ;

				auto	p { Canon::cache.try_emplace (obs, code) } ;
				if (p.second && blab > 1)
				    {
				    cout << "cache add: " << obs << " -> "
					 << (sgn1 * sgn2 < 0 ? "-" : "")
					 << *this << "\n" << flush ;
				    }
				}
			    }
			}
		    else if (blab > 2) cout << "Canon: cannot cache " << *this << "\n" ;
		    }
		else if (blab > 2) cout << "Canon: cannot cache " << *this << "\n" ;
		if (blab > 1) cout << "canon S1 returning " << sgn0 * sgn1
			       << " " << *this << "\n" << flush ;
		++global.count().cachemiss ;
		}
	    else					// cache hit
		{
		sgn1  = search->second < 0 ? -1 : 1 ;
		*this = ObsList::obs (abs(search->second)) ;
		if (blab > 1) cout << "canon S0 returning: " << sgn0 * sgn1
			       << " " << *this << "\n" << flush ;
		++global.count().cachehits ;
		}
	    }
	else						// compare starting chunks
	    {
	    int  	maxtry   { (isFf || !mid) ? 2 : 4 } ;
	    uint 	minscore { UINT_MAX } ;
	    auto&	symmset  { Canon::symmset } ;

	    if (blab > 3) cout << "\nspec obs " << *this << " mid " << mid
			       << " maxtry " << maxtry << "\n" << flush ;

	    vector<pair<int,int>> trial ;

	    for (int shot(0) ; shot < maxtry ; ++shot)
		{
		bool	reverse ( shot & 1 ) ;
		bool	flip  { isFf && reverse } ;
		int	start ( shot > 1 ? mid : (flip ? size()-1 : 0) ) ;
		uint	indx  { Canon::speccode (*this, start, reverse) } ;
		Node&	node  { Canon::spectable [indx] } ;
		if (blab > 3)
		    {
		    cout << "shot " << shot << " score " << node.score
			 << " indx " << indx << " symmset " << node.indx
			 << "\n" << flush ;
		    }
		if (node.score == UINT_MAX)
		    {
		    fatal (format("bad spectable node: obs {} start {} rev {} indx {}",
			    print(), start, reverse, indx)) ;
		    }
//		SymbStr chunk { string(specchunksize,'\0') } ;
//		if (!Canon::specchunk (indx % (spectblsize()/2), chunk))
//		    {
//		    fatal (format("specchunk failed! obs {} start {} rev {} chunk {}",
//			    print(), start, reverse, chunk.print())) ;
//		    }
//		else if (blab > 3) cout << "specchunk " << chunk << "\n" << flush ;
		if (node.score < minscore)
		    {
		    trial.resize(0) ;
		    trial.emplace_back (start, node.indx) ;
		    minscore = node.score ;
		    }
		else if (node.score == minscore)
		    {
		    trial.emplace_back (start, node.indx) ;
		    }
		}

	    if (trial.size() == 1 && symmset[trial[0].second].size() == 1)
		{
		auto		symmidx { symmset[trial[0].second][0] } ;
		const auto&	symm { Symm::list [symmidx] } ;
		int		beg  { trial[0].first } ;
		int		sgn2 { trans (symm, beg) } ;
		if (blab > 2)
		    {
		    cout << "spec best: " << symm.name
			 << " at " << beg << " -> " << sgn0 * sgn2
			 << " * " << *this << "\n" << flush ;
		    }
		sgn1 = sgn2 ;
		}
	    else				// compare remaining possibilities
		{
		int bestsymm {-1} ;
		int beststart ;
		Obs best { *this } ;

		for (auto [start,setindx] : trial)
		    {
		    for (auto symmidx : symmset[setindx])
			{
			const auto&	symm { Symm::list [symmidx] } ;
			Obs  		obs  { *this } ;
			int  		sgn2 { obs.trans (symm, start) } ;

			if (blab > 2)
			    {
			    cout << "spec try: " << symm.name
				 << " at " << start << " -> " << sgn0 * sgn2
				 << " * " << obs << "\n" << flush ;
			    }
			if (bestsymm < 0 || memcmp (best.data(), obs.data(), size()) > 0)
			    {
			    bestsymm  = symmidx ;
			    beststart = start ;
			    best = obs ;
			    sgn1 = sgn2 ;
			    }
			}
		    }
		if (blab > 2)
		    {
		    cout << "spec best: " << Symm::list[bestsymm].name
			 << " at " << beststart << " -> " << sgn0 * sgn1
			 << " * " << best << "\n" << flush ;
		    }
		*this = best ;
		}
	    if (blab > 1) cout << "canon S2 returning " << sgn0 * sgn1
			   << " " << *this << "\n" << flush ;
	    }
	return sgn0 * sgn1 ;
	}
    }

bool Canon::loopchunk (uint code, SymbStr& s)	// Map looptbl index -> loop chunk
    {
    int	 dim { theory.dim } ;
    int	 base { 2 * dim } ;
    symb x { X } ;

    for (auto p { s.rbegin() } ; p < s.rend() ; code /= base)
	{
	symb y ( code % base ) ;
	if (y >= dim) y = refl(y - dim) ;
	if (x != X && ligature(y,x)) return false ;
	*p++ = x = y ;
	}
    return true ;
    }

bool Canon::specchunk (uint code, SymbStr& s)	// Map spectbl index -> spec chunk
    {
    int	 dim  { theory.dim } ;
    int	 base { 2 * dim } ;
    int  len  ( s.size() ) ;				// must equal specchunksize
    uint blab { Blab::blablevel[BLAB::CANON] } ;
    if (blab > 3) cout << "\t specchunk: code " << code << "\n" ;

    for (auto p { s.rbegin() } ; p < s.rend() ; code /= base)	// lay down links
	{
	symb y ( code % base ) ;
	if (y >= dim) y = refl(y - dim) ;
	*p++ = y ;
	}

    if (!theory.euclid)					// add E's
	{
	int	b    { theory.nf ? 6 : 4 } ;
	short	xtra ( code % b ) ;
	short	midE ( code / b ) ;

	if (midE)
	    {
	    switch (xtra)
		{
		case 0: s[0] |= 0x08 ; s[midE] |= 0x08 ; break ;	// E_0,     E_midE
		case 1: s[0] |= 0x08 ; s[midE] |= 0x10 ; break ;	// E_0,     Elink_midE
		case 2: s[0] |= 0x10 ; s[midE] |= 0x08 ; break ;	// Elink_0, E_midE
		case 3: s[0] |= 0x10 ; s[midE] |= 0x10 ; break ;	// Elink_0, Elink_midE
		case 4: s[midE] |= 0x08 ; break ;			// link_0,  E_midE, 
		case 5: s[midE] |= 0x10 ; break ;			// link_0,  Elink_midE
		}
	    }
	else 
	    {
	    switch (xtra)
		{
		case 0: s[0] |= 0x18 ; break ;				// EE_0
		case 1: s[0] |= 0x20 ; break ;				// EElink_0
		case 2: s[0] |= 0x08 ; break ;				// E_0
		case 3: s[0] |= 0x10 ; break ;				// Elink_0
		case 4: break ;						// link_0
		case 5: return false ;
		}
	    }
	}
    if (blab > 3) cout << "\t specchunk: -> " << s << "\n" ;

    symb x { X } ;						// check ligatures
    for (auto p = s.begin() ; p < s.end() ; x = *p++)
	{
	if (ligature(x,*p) == X) fatal (format("Bad ligature {}", s.print())) ;
	if (x != X && ligature(x,*p)) return false ;
	}
    return true ;
    }

uint Canon::speccode (const Obs& obs, int start, bool reverse)	// Map spec chunk -> spectbl index
    {
    int  dim  { theory.dim } ;
    int  len  ( obs.size() ) ;
    int  inc  { reverse ? -1 : +1 } ;
    int  k    { start } ;
    uint blab { Blab::blablevel[BLAB::CANON] } ;
    if (blab > 3) cout << "speccode " << obs << " len " << len
			    << " start " << start << " reverse " << reverse << "\n" ;

    if (obs.is_fermi())						// skip fermions
	{
	len -= 2 ;
	k = reverse ? len : 1 ;
	}
    symb c    ( reverse ? conj(obs[k]) : obs[k] ) ;
    int  x    { direction(c) } ;
    uint indx ( x >= dim ? refl(x) + dim : x ) ;
    uint xtra (0) ;

    if (!theory.euclid)					// initial E or EE?
	{
	if      (isEE(c))	xtra = 0 ;
	else if (isEElink(c))	xtra = 1 ;
	else if (isE(c))	xtra = 2 ;
	else if (isElink(c))	xtra = 3 ;
	else			xtra = 4 ;
	if (blab > 3) cout << "speccode A: k " << k << " " << symbname[c]
				<< " indx " << indx << " xtra " << xtra << "\n" ;

	int b { theory.nf ? 6 : 4 } ;
	for (int m(1) ; m < specchunksize ; ++m)
	    {
	    k += inc ;
	    if (k < 0)		k += len ;
	    else if (k >= len)	k -= len ;

	    c = reverse ? conj(obs[k]) : obs[k] ;		// convert link -> code
	    x = direction(c) ;
	    indx = 2 * dim * indx + (x >= dim ? refl(x) + dim : x) ;

	    if      (isE(c))	 xtra = m * b + 2 * xtra - 4 ;	// encode mid-E
	    else if (isElink(c)) xtra = m * b + 2 * xtra - 3 ;

	    if (blab > 3) cout << "speccode B: k " << k << " " << symbname[c]
				    << " indx " << indx << " xtra " << xtra << "\n" ;
	    }
	}
    else  // theory.euclid == true
	{
	for (int m(1) ; m < specchunksize ; ++m)
	    {
	    k += inc ;
	    if (k < 0)		k += len ;
	    else if (k >= len)	k -= len ;

	    c = reverse ? conj(obs[k]) : obs[k] ;		// convert link -> code
	    x = direction(c) ;
	    indx = 2 * dim * indx + (x >= dim ? refl(x) + dim : x) ;
	    }
	}
    indx += xtra * ipow(2*dim,specchunksize) + reverse * spectblsize()/2 ;

    if (blab > 3) cout << "speccode C: indx " << indx << " xtra " << xtra
		       << " tblsize " << spectblsize() << "\n" ;
    else if (blab > 3) cout << "speccode = " << indx << "\n" << flush ;

    if (indx >= spectblsize()) fatal ("speccode: bad indx") ;
    return indx ;
    }

void Canon::looptblinit()					// Initialize looptable
    {
    int		nodecount(0) ;
    SymbStr	chk (string(loopchunksize,'\0')) ;

    for (int code(0) ; code < looptable.size() ; ++code)	// run over all loop chunks
	{
	Node&	node {looptable[code]} ;
	SymmSet	symmlist ;
	string	symmstr ;

	if (!loopchunk (code, chk)) continue ;

	for (uint symmnum(0) ; symmnum < Symm::list.size() ; ++symmnum)
	    {
	    const auto&	symm { Symm::list[symmnum] } ;
	    uint	score (0) ;

	    if (symm.isCodd())
		{
		for (auto p = chk.crbegin() ; p < chk.crend() ; ++p)
		    {
		    score = (score << loopbits) + symm.map[*p] ;
		    }
		}
	    else
		{
		for (auto p = chk.cbegin() ; p < chk.cend() ; ++p)
		    {
		    score = (score << loopbits) + symm.map[*p] ;
		    }
		}
	    if (score < node.score)
		{
		node.score = score ;
		symmlist.clear() ;
		symmlist.push_back (symmnum) ;
		symmstr = Symm::list[symmnum].name ;
		//symmstr = std::to_string(symmnum) ;
		}
	    else if (score == node.score)
		{
		symmlist.push_back (symmnum) ;
		symmstr += " " + Symm::list[symmnum].name ;
		//symmstr += " " + std::to_string(symmnum) " ;
		}
	    }
	symmstr.shrink_to_fit() ;
	symmlist.shrink_to_fit() ;
	node.indx = symmset.store (symmstr, symmlist) ;
	++nodecount ;
	}
    stats.looptblnodes = nodecount ;
    }

void Canon::spectblinit()					// Initialize spectable
    {
    int		nodecount(0) ;
    uint	halftbl	( spectable.size()/2 ) ;
    uint	csymm	( Symm::list.size()/2 ) ;
    SymbStr	chk	(string(specchunksize,'\0')) ;

    for (int code(0) ; code < halftbl ; ++code)
	{
	if (!specchunk (code, chk)) continue ;

	for (int rev(0) ; rev < 2 ; ++rev)
	    {
	    Node&	node    { spectable[code + rev * halftbl] } ;
	    uint	symmnum { rev ? csymm : 0 } ;
	    uint	symmend ( rev ? Symm::list.size() : csymm ) ;
	    SymmSet	symmlist ;
	    string	symmstr ;

	    for (; symmnum < symmend ; ++symmnum)
		{
		const auto&	symm { Symm::list[symmnum % csymm] } ;
		uint 		score (0) ;

		for (auto p = chk.cbegin() ; p < chk.cend() ; ++p)
		    {
		    score = (score << specbits) + symm.map[*p] ;
		    }
		if (score < node.score)
		    {
		    node.score = score ;
		    symmlist.clear() ;
		    symmlist.push_back (symmnum) ;
		    symmstr = Symm::list[symmnum].name ;
		    //symmstr = std::to_string(symmnum) ;
		    }
		else if (score == node.score)
		    {
		    symmlist.push_back (symmnum) ;
		    symmstr += " " + Symm::list[symmnum].name ;
		    //symmstr += " " + std::to_string(symmnum) " ;
		    }
		}
	    symmstr.shrink_to_fit() ;
	    symmlist.shrink_to_fit() ;
	    node.indx = symmset.store (symmstr, symmlist) ;
	    ++nodecount ;
	    if (Blab::blablevel[BLAB::CANON] > 2)
		{
		cout << "spectbl node[" << code + rev * halftbl << "]: "
			  << " code " << code << " chunk " << chk.print()
			  << " symmset " << node.indx << " " << symmstr << "\n" ;
		}
	    }
	}
    stats.spectblnodes = nodecount ;
    stats.symmsubsets  = symmset.size() ;
    }

void CanonCache::load (const Obsset& inbox)	// Load short inbox Obs into cache
    {
    freeze = false ;
    for (const Obs& obs : inbox)
	{
	if (obs.is_Loop()  && obs.size() > loopchunksize) continue ;
	if (!obs.is_Loop() && obs.size() > specchunksize) continue ;
	Obs o {obs} ;
	o.canon() ;
	}
    freeze = true ;
    }

void CanonCache::reload ()			// Reload cache
    {
    freeze = false ;
    for (const auto* ptr : ObsList::obs)
	{
	if (ptr->is_Loop()  && ptr->size() > loopchunksize) continue ;
	if (!ptr->is_Loop() && ptr->size() > specchunksize) continue ;
	Obs o {*ptr} ;
	o.canon() ;
	}
    freeze = true ;
    }

ostream& operator<< (ostream& stream, const CanonCache& cache)
    {
    stream << "Short Obs cache:\n" ;
    for (const auto& [key,indx] : cache)
	{
	const Obs& a { ObsList::obs (abs(indx)) } ;
	stream	<< std::left  << std::setw(12) << SymbStr(key).print() << " -> "
		<< std::right << std::setw(2) << indx << ": " << a << "\n" ;
	}
    return stream ;
    }
