#include "Gen.h"
#include "Rep.h"
#include "Global.h"
#include "Commute.h"
#include "Gripe.h"
#include "Numerics.h"
#include "Blab.h"

Gen::Gen (const Op& op)				// Construct w/o projection
    :						// (only for Test::jacobi)
    type  (op.type), T_odd (op.is_coord())
    {
    emplace_back (OpTerm (Op::store (op), 1)) ;
    inner_commute () ;
    }

Gen::Gen (const OpSum& s, GenHdr& hdr, doub coeff)	// Constructor
    :
    OpSum (s), type  ((OpType) hdr.type), order (hdr.order), T_odd (hdr.T_odd),
    imag  ((theory.euclid ? (hdr.len % 2) : hdr.T_odd) ^ !Rep::list[hdr.rep].C_even),
    coeff (coeff)
    { inner_commute() ; }

Gen::Gen (const Op& op, const Proj& proj)		// Construct with given projector
    :
    order   (op.order), T_odd   (op.is_coord()),
    imag    ((theory.euclid ? op.is_FermionO() : op.is_coord()) ^ !proj.C_even())
    {
    if (Blab::blablevel[BLAB::GEN] > 1)
	{
	cout << "Gen::Gen (Op&, Proj&) " << op << " " << proj.name << "\n" ;
	}
    for (uint i(0) ; i < proj.size() ; ++i)
	{
	if (proj[i])
	    {
	    auto [sgn,trans]	{ Symm::list[i](op) } ;
	    uint indx		{ Op::store (trans) } ;
	    emplace_back (OpTerm {indx, proj[i] * sgn} ) ;
	    }
	}
    if (int k = OpSum::collect(true))		// collect terms
	{
	coeff = k/proj.denom ;
	type = op.type ;
	inner_commute() ;
	}
    }

Gen::Gen (const OpSum& s, const Proj& proj)		// Construct with given projector
    {
    if (Blab::blablevel[BLAB::GEN] > 1)
	{
	cout << "Gen::Gen (OpSum&, proj&) " << proj.name << "\n" ;
	}
    OpType	tmptype { OpType::Invalid } ;
    bool	iseuc	{ theory.euclid } ;
    bool	Todd ;

    for (const OpTerm& t : s)
	{
	const Op op { Op::list[t.item] } ;

	if (tmptype == OpType::Invalid)
	    {
	    tmptype = op.type ;
	    order   = op.order ;
	    T_odd   = op.is_coord() ;
	    imag    = (iseuc ? op.is_FermionO() : op.is_coord()) ^ !proj.C_even() ;
	    }
	else if (tmptype != op.type  || T_odd  != op.is_coord() || order != op.order)
	    {
	    cout << "Gen::Gen: gen " << *this << " + op " << op << "\n" ;
	    if (op.type  != tmptype)	gripe ("Differing OpTypes") ;
	    if (op.order != order)	gripe ("Differing Op orders") ;
	    if (op.is_coord() != T_odd) gripe ("Differing T-symms") ;
	    }

	for (uint i(0) ; i < proj.size() ; ++i)
	    {
	    if (proj[i])
		{
		auto [sgn,trans] { Symm::list[i](op) } ;
		uint indx	 { Op::store (trans) } ;
		emplace_back (indx, t.coeff * proj[i] * sgn) ;
		}
	    }
	}
    if (int k = OpSum::collect(true))		// collect terms
	{
	coeff = k/proj.denom ;
	type  = tmptype ;
	inner_commute() ;
	}
    }

void Gen::settype (Op& op)			// Set & check Gen type
    {
    if (Blab::blablevel[BLAB::GEN] > 1)
	{
	cout << "Gen::settype gen " << *this << " Op " << op
	     << " size " << size() << "\n" ;
	}
    if (type == OpType::Invalid)
	{
	type  = op.type ;
	T_odd = op.is_coord() ;
	}
    else if (op.type != type || op.is_coord() != T_odd)
	{
	cout << "Gen::settype: gen " << *this << " + op " << op << "\n" ;
	if (op.type       != type)	gripe ("Differing OpTypes") ;
	if (op.is_coord() != T_odd)	gripe ("Differing T-symms") ;
	}
    }

int Gen::project (Op& op)			// Project onto all reps
    {
    int added(0) ;
    for (int repnum(0) ; repnum < Rep::list.size() ; ++repnum)
	{
	for (const Proj& proj : Rep::list[repnum])
	    {
	    Gen  tmp { op, proj } ;

	    if (tmp.valid() && isnew (repnum, tmp))
		{
		auto& info  { global.info(tmp.is_Fermion()) } ;
		auto& gens  { info.gens[repnum] } ;
		auto& neven { info.neven[repnum] } ;
		int   indx  ( tmp.T_odd ? gens.size() : neven ) ;

		if (Gen::gennorm) tmp.normalize (repnum) ;
		if (tmp.T_odd) gens.push_back (tmp) ;
		else gens.insert (gens.begin() + neven++, tmp) ;

		if (tmp.order > global.maxgen()) global.maxgen() = tmp.order ;
		++added ;

		if (Blab::blablevel[BLAB::GEN])
		    {
		    cout << proj.name << " gen " << global.fg()
			 << indx << (tmp.T_odd ? "*" : "")
			 << " (" << tmp.order << ") " << tmp
			 << "\n\treduction: " << tmp.reduction << "\n";
		    }
		}
	    }
	}
    return added ;
    }

int Gen::project (OpSum&& s)			// Project onto all reps
    {
    int added(0) ;
    for (int repnum(0) ; repnum < Rep::list.size() ; ++repnum)
	{
	for (const auto& proj : Rep::list[repnum])
	    {
	    Gen	tmp { s, proj } ;

	    if (tmp.valid() && isnew (repnum,tmp))
		{
		auto& info  { global.info(tmp.is_Fermion()) } ;
		auto& gens  { info.gens[repnum] } ;
		auto& neven { info.neven[repnum] } ;
		int   indx  ( tmp.T_odd ? gens.size() : neven ) ;

		if (Gen::gennorm) tmp.normalize (repnum) ;
		if (tmp.T_odd) gens.push_back (tmp) ;
		else gens.insert (gens.begin() + neven++, tmp) ;

		if (tmp.order > global.maxgen()) global.maxgen() = tmp.order ;
		++added ;

		if (Blab::blablevel[BLAB::GEN] > 1)
		    {
		    cout << proj.name << " gen #" << global.fg()
			 << indx << (tmp.T_odd ? "*" : "")
			 << " (" << tmp.order << ") " << tmp
			 << "\n\treduction: " << tmp.reduction << "\n";
		    }
		}
	    }
	}
    return added ;
    }

bool Gen::isnew (int repnum, const Gen& b)		// New generator?
    {
    bool	stage { b.is_Fermion() } ;
    const auto& gens  { global.info(stage).gens[repnum] } ;
    uint	ngen  ( gens.size() ) ;
    uint	start { stage ? global.nopG() : 0 } ;
    int		maxop (-1) ;

    for (const auto& a : gens)
	{
	for (const auto& s : a) if (maxop < 0 || s.item > maxop) maxop = s.item ;
	}
    for (const auto& s : b)
	{
	if (maxop < 0 || s.item > maxop) return true ;
	}

    Nummtx genmtx (maxop-start+1, ngen+1) ;

    for (int i(0) ; i < ngen ; ++i)
	{
	for (const auto& s : gens[i]) genmtx(s.item - start,i) = s.coeff ;
	}
    for (const auto& s : b) genmtx(s.item - start,ngen) = s.coeff ;

    return arma::rank (genmtx) == genmtx.n_cols ;
    }

bool Gen::allzero () const				// Vanishing Gen?
    {
    for (const auto t : *this) if (t.coeff) return false ;
    return true ;
    }

Gen& Gen::collect ()					// Collect terms
    {
    if (int k = OpSum::collect()) coeff *= k ;
    if (allzero()) resize (0) ;
    return *this ;
    }

void Gen::normalize (int repnum)			// Normalize generator
    {
    const auto&		Hterms	{ global.info().Hterms } ;
    ObsList		tmplist	{ "GenNorm" } ;
    doub		normsq	(0) ;

    for (ushort i(0) ; i < Hterms.size() ; ++i)
	{
	const ObsPoly&	poly { repnum ? Hterms[i].poly : Hterms[i].cpoly } ;
	const ObsList&	list { poly.obslist() } ;

	for (const auto& term : poly)
	    {
	    const auto& obs { list(term[0]) } ;
	    if (obs.is_EEloop()  && obs.size() == 1 ||
		obs.is_Entropy() && !is_Fermion()   ||
		obs.is_Fermion() && obs.size() == 2 && theory.euclid)
		{
		PolyMap	tmp { tmplist } ;
		PolyMap	ans { tmplist } ;
		Commute::commute_poly (*this, poly, tmp) ;
		Commute::commute_poly (*this, tmp,  ans) ;
		if (imag) ans.negate() ;
		if (theory.euclid && is_Fermion())
		    {
		    PolyTerm reindx { tmplist.catalog (obs) } ;
		    normsq += abs(ans[reindx]) ;
		    }
		else normsq += ans[Polyindx()] ;
		}
	    }
	}
    if (!normsq)
	{
	for (const auto& opterm : *this)
	    {
	    doub x { coeff * opterm.coeff } ;
	    normsq += x * x ;
	    }
	if (Blab::blablevel[BLAB::GEN])
	    cout << "Gen::normalize: method 2 for " << *this << "\n" ;
	}
    coeff /= sqrt (normsq) ;
    }

void Gen::normalize ()					// Generator normalization
    {
    for (int repnum(0) ; repnum < Rep::list.size() ; ++repnum)
	{
	auto& gens { global.info().gens[repnum] } ;

	for (auto& gen : gens) gen.normalize (repnum) ;
	}
    }

void Gen::inner_commute ()			 	// Generator reduction
    {
    if (type == OpType::Loop) return ;

    uint blab { Blab::blablevel[BLAB::GEN] } ;
    if (blab > 1) cout << "--- inner_commute: " << *this << " " << (int)type << "\n" ;

    PolyMap map { reduction.obslist() } ;

    if (type == OpType::Eloop)
	{
	for (const auto& opterm : *this)
	    {
	    Obs		eloop  { Op::list[opterm.item] } ;
	    PolyTerm	factor { Polyindx(), opterm.coeff } ;
	    Commute::do_inner (eloop, factor, reduction.obslist(), map) ;
	    }
	}
    else if (type == OpType::Fermion)
	{
	ObsList& list { reduction.obslist() } ;
	for (const auto& opterm : *this)
	    {
	    Op& op { Op::list[opterm.item] } ;
	    if (!op.staggered() || !op.isclosed()) continue ;
	    int sgn  { isstag(op.front()) ? 1 : -1 } ;
	    Obs loop { op.begin()+1, op.end()-1, ObsType::Loop } ;
	    PolyTerm result { list.assess(loop) * opterm.coeff * sgn } ;
	    if (result.coeff) map.add (std::move(result)) ;
	    }
	}
    reduction.clear() ;
    reduction.push_map (map) ;

    if (blab > 2)	cout << "inner_commute: returning " << reduction << "\n" ;
    else if (blab > 1)	cout << "inner_commute: returning \n" ;
    }

int Gen::addgen (OpSum&& s)			// Add new generator
    {
    int added(0) ;
    switch (Op::list[s.front().item].type)
	{
	case OpType::Loop:
	    if (autoEgens)	added += project (Op::loop_dt(s)) ;
	    if (!theory.euclid)	added += project (std::move(s)) ; break ;
	case OpType::Eloop:
	case OpType::Fermion:	added += project (std::move(s)) ; break ;
	default:		break ;
	}
    if (added) Op::setprimary () ;
    return added ;
    }

void Gen::geninit ()				// Generator initialization
    {
    bool isham { !theory.euclid } ;
    int stage  { global.stage } ;
    char l[4]  { 'x', 'y', 'z', 'w' } ;
    char L[4]  { 'X', 'Y', 'Z', 'W' } ;
    char f[4]  { 'f', 'g', 'h', 'i' } ;
    char F[4]  { 'F', 'G', 'H', 'I' } ;

    if (stage == 0)
	{
	OpType loop { OpType::Loop } ;

	for (int i(0) ; i < theory.dim ; ++i)		// 1x1 plaq
	    {
	    for (int j(i) ; ++j < theory.dim ;)
		{
		Op plaq { string {l[i],l[j],L[i],L[j]}, loop, 2 } ;
		if (autoEgens)	project (Op::loop_dt (plaq)) ; 
		if (isham)	project (plaq) ;
		}
	    }

	for (int i(0) ; i < theory.dim ; ++i)		// 2x1 rectangle, possibly bent
	    {
	    for (int j(0) ;  j < theory.dim ; ++j)
		{
		for (int k(0) ; k < theory.dim ; ++k)
		    {
		    if (i == k || j == k) continue ;
			{
			Op rect { string {l[i],l[j],l[k],L[j],L[i],L[k]}, loop, 4 } ;
			if (autoEgens)	project (Op::loop_dt (rect)) ; 
			if (isham)	project (rect) ;
			}
		    }
		}
	    }

	for (int i(0) ; i < theory.dim ; ++i)		// 2x1 fig-8, possibly bent
	    {
	    for (int j(0) ;  j < theory.dim ; ++j)
		{
		for (int k(0) ;  k < theory.dim ; ++k)
		    {
		    if (i == k || j == k) continue ;
		    Op fig8 { string {l[i],l[k],l[j],L[k],L[j],l[k],L[i],L[k]}, loop, 4 } ;
		    if (autoEgens)	project (Op::loop_dt (fig8)) ; 
		    if (isham)		project (fig8) ;
		    }
		}
	    }

	for (int i(0) ; i < theory.dim ; ++i)		// 1x1 plaq^2
	    {
	    for (int j(i) ;  j < theory.dim ; ++j)
		{
		if (i == j) continue ;
		Op plaq2 { string {l[i],l[j],L[i],L[j],l[i],l[j],L[i],L[j]}, loop, 4 } ;
		if (autoEgens)	project (Op::loop_dt (plaq2)) ; 
		if (isham)	project (plaq2) ;
		}
	    }

	for (int i(0) ; i < theory.dim ; ++i)		// xx..x polyakov
	    {
	    if (theory.box.comp[i])
		{
		Op polyakov  { string (theory.box.comp[i],   l[i]), loop, 2 } ;
		if (autoEgens)	project (Op::loop_dt (polyakov)) ; 
		if (isham)	project (polyakov) ;
		}
	    }

	for (int i(0) ; i < theory.dim ; ++i)		// plaq.polyakov
	    {
	    if (theory.box.comp[i])
		{
		for (int j(0) ; j < theory.dim ; ++j)
		    {
		    if (i == j) continue ;
		    string plaq {l[i],l[j],L[i],L[j]} ;
		    string poly (theory.box.comp[i], l[i]) ;
		    Op plaqpoly { plaq+poly, loop, 4 } ;
		    if (autoEgens)	project (Op::loop_dt (plaqpoly)) ; 
		    if (isham)		project (plaqpoly) ;
		    }
		}
	    }
	}
    else if (theory.nf > 0)
	{
	OpType ferm { OpType::Fermion } ;

	for (int k(0) ; k < theory.nf ; k += 2)		// Ff, Gf
	    {
	    if (!isham)
		{
		Op Gf { string {F[k+1],f[k]}, ferm, 0 } ;
		project (Gf) ;
		}
	    }

	for (int k(0) ; k < theory.nf ; k += 2)		// Fxf, Gxf
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		Op Gxf { string {F[k+1],l[i],f[k]}, ferm, 1 } ;
		project (Gxf) ;
		if (!isham) continue ;
		Op Fxf { string {F[k],l[i],f[k]}, ferm, 1 } ;
		project (Fxf) ;
		}

	for (int k(0) ; k < theory.nf ; k += 2)		// Fxyf, Gxyf
	    for (int i(0) ; i < theory.dim ; ++i)
		for (int j(0) ;  j < theory.dim ; ++j)
		    {
		    Op Gxyf { string {F[k+1],l[i],l[j],f[k]}, ferm, 2 } ;
		    project (Gxyf) ;
		    if (!isham) continue ;
		    Op Fxyf { string {F[k],l[i],l[j],f[k]}, ferm, 2 } ;
		    project (Fxyf) ;
		    }
	
	for (int k(0) ; k < theory.nf ; k += 2)		// FxyXYf, GxyXYf
	    for (int i(0) ; i < theory.dim ; ++i)
		for (int j(i) ; ++j < theory.dim ;)
		    {
		    Op Gplaqf { string {F[k+1],l[i],l[j],L[i],L[j],f[k]}, ferm, 4 } ;
		    project (Gplaqf) ;
		    if (!isham) continue ;
		    Op Fplaqf { string {F[k],l[i],l[j],L[i],L[j],f[k]}, ferm, 4 } ;
		    project (Fplaqf) ;
		    }
	}
    }

ostream& Gen::print (ostream& stream) const		// Short form print
    {
    int k(0) ;
    coeffprt (stream, coeff) ;
    stream << (imag ? "i " : "") << "( " ;
    for (auto t : *this)
	{
	if (t.coeff != 0)
	    {
	    ++k ;
	    coeffprt (stream, t.coeff); 
	    stream << Op::list[t.item] << " + ..." ;
	    break ;
	    }
	}
    stream << (k ? " )" : " 0 )") ;
    return stream ;
    }

ostream& operator<< (ostream& stream, const Gen& gen)
    {
    int k(0) ;
    coeffprt (stream, gen.coeff) ;
    stream << (gen.imag ? "i " : "") << "(" ;
    for (auto t : gen)
	{
	if (t.coeff != 0)
	    {
	    stream << "\n\t" ;
	    coeffprt (stream, t.coeff); 
	    stream << Op::list[t.item] ;
	    ++k ;
	    }
	}
    if (!k) stream << " 0" ;
    stream << " )" ;
    if (gen.reduction.size())
	{
	stream << "\n reduction = " ;
	coeffprt (stream, gen.coeff) ;
	stream << "(" <<  gen.reduction << ")" ;
	}
    return stream ;
    }
