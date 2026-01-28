#include "Build.h"
#include "Theory.h"
#include "Global.h"
#include "Poly.h"
#include "Canon.h"
#include "Commute.h"
#include "Numerics.h"
#include "Rep.h"
#include "Save.h"
#include "Blab.h"
#include "Gripe.h"

static int		cord ;			// Current creation order
static Obsset		newobs ;		// Newly generated Obs
static std::mutex	obsmutex ;		// newobs insertion lock
constexpr int		single_thread = 1024 ;	// multi-threading threshold

void Build::mk_obs (int target)			// Build observables
    {
    ObsList::freeze = false ;
    ObsList::obs.approx = false ;
    Canon::cache.freeze = false ;
    if (global.stage == Global::Gauge)
	{
	auto& maxord { global.maxord() } ;
	while (maxord < target)
	    {
	    if (++maxord % 4) continue ;
	    mk_loops () ;
	    if (!theory.euclid) { mk_Eloops() ; mk_EEloops() ; }
	    }
	}
    else // global.stage == Global::Fermi
	{
	auto& maxord { global.maxord() } ;
	while (maxord < target)
	    {
	    if (++maxord % 2) continue ;
	    mk_fermions () ;
	    if (!theory.euclid) mk_Efermions() ;
	    }
	ObsList::obs.do_fermi_init () ;
	}
    Obsset().swap(newobs) ;
    Obsset().swap(ObsList::inbox) ;
    ObsList::freeze = true ;
    ObsList::obs.shrink_to_fit() ;
    ObsList::obs.approx = global.approx ;
    Canon::cache.freeze = true ;
    Global::mk_bcktlist (global.stage) ;
    Numerics::numericsinit () ;

    if (Blab::blablevel[BLAB::BUILD])
	cout << "Total # Obs:\t" << global.nobs() << "\n" << flush ;
    }

void Build::mk_loops ()				// Build Loop's
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	bckt	{ global.info().bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
    uint	numobs	{ obslist.nobs() } ;

    for (cord = maxord/2 ; cord <= maxord - 4 ; cord += 2)
	{
	if (blab > 2) cout << "\nmk_loops: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.nobs() ) ;
	Global::mk_bcktlist (Global::Gauge) ;
	if (global.nobs() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_Loop_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.nobs() + newobs.size() ;
	global.nobsG() += delta_n ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    Obsset().swap(newobs) ;
	    }
	if (blab > 1 && delta_n)
	    cout << "\n" << std::left << std::setw(12) << "Loops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "#obs = " << obslist.nobs() << "\n" << flush ;
	}
    if (blab == 1 && numobs < obslist.nobs())
	cout << "Loops     (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.nobs() - numobs << " added, #obs = "
	     << obslist.nobs() << "\n" << flush ;
    }

void Build::mk_Eloops ()			// Build Eloop's
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	bckt	{ global.info().bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
    uint	numobs	{ obslist.nobs() } ;
    int		maxcord	{ maxord > 4 ? maxord - 4 : 2 } ;

    for (cord = maxord/2 ; cord <= maxcord ; cord += 2)
	{
	if (blab > 2) cout << "\nmk_Eloops: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.nobs() ) ;
	Global::mk_bcktlist (Global::Gauge) ;
	if (global.nobs() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_Eloop_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.nobs() + newobs.size() ;
	global.nobsG() += delta_n ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    Obsset().swap(newobs) ;
	    }
	if (blab > 1 && delta_n)
	    cout << "\n" << std::left << std::setw(12) << "Eloops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "#obs = " << obslist.nobs() << "\n" << flush ;
	}
    if (blab == 1 && numobs < obslist.nobs())
	cout << "Eloops    (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.nobs() - numobs << " added, #obs = "
	     << obslist.nobs() << "\n" << flush ;
    }

void Build::mk_EEloops ()			// Build EEloop's
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	bckt	{ global.info().bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
    uint	numobs	{ obslist.nobs() } ;

    for (cord = maxord/2 - 2 ; cord <= maxord - 4 ; cord += 2)
	{
	if (blab > 2) cout << "\nmk_EEloops: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.nobs() ) ;
	Global::mk_bcktlist (Global::Gauge) ;
	if (global.nobs() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_EEloop_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.nobs() + newobs.size() ;
	global.nobsG() += delta_n ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    Obsset().swap(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "EEloops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "#obs = " << obslist.nobs() << "\n" << flush ;
	}
    if (blab == 1 && numobs < obslist.nobs())
	cout << "EEloops   (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.nobs() - numobs << " added, #obs = "
	     << obslist.nobs() << "\n" << flush ;
    }

void Build::mk_fermions ()			// Build Fermion's
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	bckt	{ global.info().bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
    uint	numobs	{ obslist.nobs() } ;

    for (cord = maxord/2 ; cord <= maxord - 2 ; ++cord)
	{
	if (blab > 2) cout << "mk_fermions: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.nobs() ) ;
	Global::mk_bcktlist (Global::Fermi) ;
	if (global.nobs() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_Fermion_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.nobs() + newobs.size() ;
	global.nobsF() += delta_n ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    Obsset().swap(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "Fermions"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "#obs = " << obslist.nobs() << "\n" << flush ;
	}
    if (blab == 1 && numobs < obslist.nobs())
	cout << "Fermions  (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.nobs() - numobs << " added, #obs = "
	     << obslist.nobs() << "\n" << flush ;
    }

void Build::mk_Efermions ()			// Build Efermion's
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	bckt	{ global.info().bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
    uint	numobs	{ obslist.nobs() } ;

    for (cord = maxord/2 ; cord <= maxord - 2 ; ++cord)
	{
	if (blab > 2) cout << "mk_Efermions: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.nobs() ) ;
	Global::mk_bcktlist (Global::Fermi) ;
	if (global.nobs() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_Efermion_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.nobs() + newobs.size() ;
	global.nobsF() += delta_n ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    Obsset().swap(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "Efermions"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "#obs = " << obslist.nobs() << "\n" << flush ;
	}
    if (blab == 1 && numobs < obslist.nobs())
	cout << "Efermions (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.nobs() - numobs << " added, #obs = "
	     << obslist.nobs() << "\n" << flush ;
    }

void Build::mk_grad()				// Build gradient
    {
    if (global.interrupt) return ;

    auto&	info	{ global.info() } ;
    auto&	grad	{ global.data().grad } ;
    const auto&	gens	{ info.gens.front() } ;
    const auto&	Hterms	{ info.Hterms } ;
    ushort	ngens	{ info.neven.front() } ;
    ushort	nterms	( Hterms.size() ) ;
    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    PolyMap	ans	{ ObsList::obs } ;

    if (!global.maxord()) gripe ("Make some observables first!") ;

    grad.clear() ;
    if (blab) cout << "Gradient:    " << flush ;
    for (ushort i(0) ; i < nterms ; ++i)
	{
	const ObsPoly& poly { Hterms[i].cpoly } ;

	for (ushort j(0) ; j < ngens ; ++j)
	    {
	    if (!gens[j].T_odd)
		{
		Commute::commute_poly (gens[j], poly, ans) ;
		if (Hterms[i].imag && gens[j].imag) ans.negate() ;
		}
	    grad.add (Poly {GradHdr {j,i}}, ans.negate()) ;
	    ans.clear() ;
	    }
	}
    if (grad.entry().id != RecordID::Grad) fatal ("mk_grad: bad record ID!") ;
    grad.shrink_to_fit() ;
    grad.entry().ncol   = nterms ;
    grad.entry().nrow   = ngens ;
    grad.entry().reclen = grad.size() ;
    if (blab) cout << "\t\tdone\n" << flush ;
    }

void Build::mk_curv (int repnum)			// Build curvature
    {
    if (global.interrupt) return ;

    auto&	info	{ global.info() } ;
    auto&	curv	{ global.data().curv[repnum] } ;
    const auto& repnam	{ Rep::list[repnum].name } ;
    const auto&	gens	{ info.gens[repnum] } ;
    const auto&	Hterms	{ info.Hterms } ;
    ushort	nterms	( Hterms.size() ) ;
    ushort	ngens   ( gens.size() ) ;
    uint	blab    { Blab::blablevel[BLAB::BUILD] } ;
    ObsList	tmplist { "CurvTemp" } ;
    PolyMap	tmp	{ tmplist } ;
    PolyMap	ans	{ ObsList::obs } ;

    if (!global.maxord()) gripe ("Make some observables first!") ;

    curv.clear() ;

    if (blab && ngens && nterms) cout << repnam << " curvature:    " << flush ;

    for (ushort i(0) ; i < nterms ; ++i)
	{
	const ObsPoly&	poly  { repnum ? Hterms[i].poly : Hterms[i].cpoly } ;

	for (ushort j(0) ; j < ngens ; ++j)
	    {
	    for (ushort k(0) ; k < ngens ; ++k)
		{
		if (gens[j].T_odd == gens[k].T_odd)
		    {
		    Commute::commute_poly (gens[j], poly, tmp) ;
		    Commute::commute_poly (gens[k], tmp,  ans) ;
		    if (gens[j].imag != gens[k].imag && !theory.euclid)
			fatal ("mk_curv: imaginary!") ;
		    else if (gens[j].imag && gens[k].imag   ||
			     gens[j].imag && Hterms[i].imag ||
			     gens[k].imag && Hterms[i].imag) ans.negate() ;
		    }
		curv.add (Poly {CurvHdr {k,j,i}}, ans) ;
		ans.clear() ;
		tmp.clear() ;
		}
	    }
	}
    if (curv.entry().id != RecordID::Curv) fatal ("mk_curv: bad record ID!") ;
    curv.shrink_to_fit() ;
    curv.entry().nslice  = nterms ;
    curv.entry().ncol    = ngens ;
    curv.entry().nrow    = ngens ;
    curv.entry().reclen  = curv.size() ;
    if (blab && ngens && nterms) cout << "\tdone\n" << flush ;
    }

void Build::mk_lagr (int repnum)			// Build Lagrange bracket
    {
    if (global.interrupt) return ;

    const auto&	gens	{ global.info().gens[repnum] } ;
    auto&	lagr	{ global.data().lagr[repnum] } ;
    const auto& repnam	{ Rep::list[repnum].name } ;
    ushort	ngens	( gens.size() ) ;
    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    PolyMap	ans	{ ObsList::obs } ;
    short	trunc	{ SHRT_MAX } ;

    if (!global.maxord()) gripe ("Make some observables first!") ;

    if (ngens && blab) cout << repnam << " Lagrange brkt: " << flush ;

    lagr.clear();
    for (ushort j(0) ; j < ngens ; ++j)
	{
	for (ushort k(0) ; k < ngens ; ++k)
	    {
	    if (gens[j].T_odd != gens[k].T_odd)
		{
		Gen newgen ;
		Commute::commute_gen (gens[j], gens[k], newgen) ;
		if (!ans.add_gen (newgen))
		    trunc = std::min(trunc, newgen.order) ;
		}
	    lagr.add (Poly {LagrHdr {k,j}}, ans) ;
	    ans.clear() ;
	    }
	}
    if (lagr.entry().id != RecordID::Lagr) fatal ("mk_lagr: bad record ID!") ;
    lagr.shrink_to_fit() ;
    lagr.entry().ncol   = ngens ;
    lagr.entry().nrow   = ngens ;
    lagr.entry().reclen = lagr.size();
    if (blab && ngens)
	{
	if (trunc == SHRT_MAX) cout << "\tdone\n" ;
	else cout << format("\ttruncation in order {} elements\n", trunc) ;
	}
    }

void Build::mk_geos()				// Build geodesic equations
    {
    if (global.interrupt) return ;

    const auto&	bckt	{ global.info().bckt } ;
    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;

    if (blab && blab < 3) cout << "geodesics:\t" << flush;
    if (!global.maxord()) gripe ("Make some observables first!") ;

    ObsList::freeze = true ;
    Canon::cache.freeze = true ;
    global.count().cleargeostats() ;
    if (global.autosave) Save::save_sys() ;
    TASK_ARENA (global.maxthread, bckt,
	FOR_EACH (bckt.begin(), bckt.end(), do_geo_bckt)) ;

    if (blab > 1) cout << "." << flush ;
    }

void Build::do_Loop_bckt (const uint3& bckt)		// Do loop build bucket 
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    ObsList&	obslist	{ ObsList::obs } ;
    PolyTerm 	zero	{ Polyindx(), 0 } ;
    PolyMap 	tmp	{ obslist } ;
    uint	numobs	{ obslist.nobs() } ;
    uint	bcktnum { bckt[0] } ;
    uint 	first	{ bckt[1] } ;
    uint 	last	{ bckt[2] } ;

    ObsList::freeze = false ;
    for (uint i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	if (obslist(i).type != ObsType::Loop)	continue ;

	Obs a	{ obslist(i) } ;
	int ord	{ a.order() } ;
	for (const auto& op : Op::list)
	    {
	    if (op.primary == false)		continue ;
	    if (op.type != OpType::Eloop)	continue ;
	    if (op.order + a.corder != cord)	continue ;
	    if (ord < maxord - 2 * op.order)	continue ;
	    if (blab > 3)
		{
		cout << "do_Loop_bckt " << bcktnum << " [" << op << "," << a
		     << "] (" << op.order << "+" << a.corder << ")\n" << flush ;
		}
	    Commute::do_commute (op, a, zero, obslist, tmp) ;
	    }
	}
    ObsList::freeze = true ;
    if (inbox.size())
	{
	std::lock_guard<std::mutex> lock (obsmutex) ;
	newobs.insert (inbox.begin(), inbox.end()) ;
	}
    if (blab > 1)
	{
	char c { (inbox.size() || obslist.size() > numobs) ? '+' : '.' } ;
	cout << c << flush ;
	}
    inbox.clear() ;
    }

void Build::do_Eloop_bckt (const uint3& bckt)		// Do Eloop build bucket
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    ObsList&	obslist	{ ObsList::obs } ;
    PolyTerm	zero	{ Polyindx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    uint	numobs	{ obslist.nobs() } ;
    uint	bcktnum { bckt[0] } ;
    uint 	first	{ bckt[1] } ;
    uint 	last	{ bckt[2] } ;

    ObsList::freeze = false ;
    for (uint i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	auto	insiz { inbox.size() } ;
	if (obslist(i).is_Loop())
	    {
	    if (obslist(i).corder != cord) continue ;

	    Op op { obslist(i) } ;
	    for (const Obs* EE : ObsList::base)
		{
		if (EE->size() != 1 || !isEE(EE->front())) continue ;
		if (blab > 3)
		    {
		    cout << "do_Eloop_bckt " << bcktnum << " inbox.size "
			 << inbox.size() << " [" << op << "," << *EE
			 << "] (" << op.order << "+" << EE->corder << ")\n" << flush ;
		    }
		Commute::do_commute (op, *EE, zero, obslist, tmp) ;
		if (blab > 3)
		    {
		    cout << "do_Eloop_bckt " << bcktnum << " inbox.size "
			 << inbox.size() << "\n" << flush ;
		    }
		}
	    }
	else if (obslist(i).is_Eloop())
	    {
	    Obs	a	{ obslist(i) } ;
	    int	ord	{ a.order() } ;
	    for (const auto& op : Op::list)
		{
		if (op.primary == false)		continue ;
		if (op.type != OpType::Eloop)		continue ;
		if (op.order + a.corder != cord)	continue ;
		if (ord < maxord - 2 * op.order)	continue ;
		if (blab > 3)
		    {
		    cout << "do_Eloop_bckt " << bcktnum << " inbox.size "
			 << inbox.size() << " [" << op << "," << a
			 << "] (" << op.order << "+" << a.corder << ")\n" << flush ;
		    }
		Commute::do_commute (op, a, zero, obslist, tmp) ;
		if (blab > 3)
		    {
		    cout << "do_Eloop_bckt " << bcktnum << " inbox.size "
			 << inbox.size() << "\n" << flush ;
		    }
		}
	    }
	}
    ObsList::freeze = true ;
    if (inbox.size())
	{
	std::lock_guard<std::mutex> lock (obsmutex) ;
	newobs.insert (inbox.begin(), inbox.end()) ;
	}
    if (blab > 1)
	{
	char c { (inbox.size() || obslist.size() > numobs) ? '+' : '.' } ;
	cout << c << flush ;
	}
    inbox.clear() ;
    }

void Build::do_EEloop_bckt (const uint3& bckt)		// Do EEloop build bckt
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    ObsList&	obslist	{ ObsList::obs } ;
    uint	numobs	{ obslist.nobs() } ;
    PolyTerm	zero	{ Polyindx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    uint	bcktnum { bckt[0] } ;
    uint 	first	{ bckt[1] } ;
    uint 	last	{ bckt[2] } ;

    ObsList::freeze = false ;
    for (uint i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	if (obslist(i).type != ObsType::EEloop) continue ;

	Obs a	{ obslist(i) } ;
	int ord	{ a.order() } ;
	for (const auto& op : Op::list)
	    {
	    if (op.primary == false)		continue ;
	    if (op.type != OpType::Eloop)	continue ;
	    if (op.order + a.corder != cord)	continue ;
	    if (ord < maxord - 2 * op.order)	continue ;
	    if (blab > 3)
		{
		cout << "do_EEloop_bckt " << bcktnum << " [" << op << "," << a
		     << "] (" << op.order << "+" << a.corder << ")\n" << flush ;
		}
	    Commute::do_commute (op, a, zero, obslist, tmp) ;
	    }
	}
    ObsList::freeze = true ;
    if (inbox.size())
	{
	std::lock_guard<std::mutex> lock (obsmutex) ;
	newobs.insert (inbox.begin(), inbox.end()) ;
	}
    if (blab > 1)
	{
	char c { (inbox.size() || obslist.size() > numobs) ? '+' : '.' } ;
	cout << c << flush ;
	}
    inbox.clear() ;
    }

void Build::do_Fermion_bckt (const uint3& bckt)		// Do Fermion build bckt
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    ObsList&	obslist	{ ObsList::obs } ;
    uint	numobs	{ obslist.nobs() } ;
    PolyTerm	zero	{ Polyindx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    uint	bcktnum { bckt[0] } ;
    uint	first	{ bckt[1] } ;
    uint	last	{ bckt[2] } ;

    ObsList::freeze = false ;
    for (uint i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	if (obslist(i).type != ObsType::Fermion) continue ;

	Obs a { obslist(i) } ;
	int ord	{ a.order() } ;
	for (const auto& op : Op::list)
	    {
	    if (op.is_coord() || !op.primary)	continue ;
	    if (op.order + a.corder != cord)	continue ;
	    if (ord < maxord - 2 * op.order)	continue ;
	    if (blab > 3)
		{
		cout << "do_Fermion_bckt " << bcktnum << " [" << op << "," << a
		     << "] (" << op.order << "+" << a.corder << ")\n" << flush ;
		}
	    Commute::do_commute (op, a, zero, obslist, tmp) ;
	    }
	}
    ObsList::freeze = true ;
    if (inbox.size())
	{
	std::lock_guard<std::mutex> lock (obsmutex) ;
	newobs.insert (inbox.begin(), inbox.end()) ;
	}
    if (blab > 1)
	{
	char c { (inbox.size() || obslist.size() > numobs) ? '+' : '.' } ;
	cout << c << flush ;
	}
    inbox.clear() ;
    }

void Build::do_Efermion_bckt (const uint3& bckt)	// Do Efermion build bckt
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    ObsList&	obslist	{ ObsList::obs } ;
    uint	numobs	{ obslist.nobs() } ;
    PolyTerm	zero	{ Polyindx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    uint	bcktnum { bckt[0] } ;
    uint 	first	{ bckt[1] } ;
    uint 	last	{ bckt[2] } ;

    ObsList::freeze = false ;
    for (uint i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	if (obslist(i).is_Fermion())
	    {
	    if (obslist(i).corder != cord) continue ;

	    Op op { obslist(i) } ;
	    if (!theory.euclid) op.front() = stag(op.front()) ;
	    for (const Obs* EE : ObsList::base)
		{
		if (EE->size() != 1 || !isEE(EE->front())) continue ;
		if (blab > 3)
		    {
		    cout << "\nmk_Efermions " << bcktnum << " [" << op << "," << *EE
			 << "] (" << op.order << "+" << EE->corder << ")\n" << flush ;
		    }
		Commute::do_commute (op, *EE, zero, obslist, tmp) ;
		}
	    }
	else if (obslist(i).is_Efermion())
	    {
	    Obs a	{ obslist(i) } ;
	    int ord	{ a.order() } ;
	    for (const auto& op : Op::list)
		{
		if (op.is_coord() || !op.primary)	continue ;
		if (op.order + a.corder != cord)	continue ;
		if (ord < maxord - 2 * op.order)	continue ;
		if (blab > 3)
		    {
		    cout << "\nmk_Efermions [" << op << "," << a << "] ("
			 << op.order << "+" << a.corder << ")\n" << flush ;
		    }
		if (op.is_Eloop())
		    Commute::do_commute  (op, a, zero, obslist, tmp) ;
		else
		    Commute::do_commuteE (op, a, zero, obslist, tmp) ;
		}
	    }
	}
    ObsList::freeze = true ;
    if (inbox.size())
	{
	std::lock_guard<std::mutex> lock (obsmutex) ;
	newobs.insert (inbox.begin(), inbox.end()) ;
	}
    if (blab > 1)
	{
	char c { (inbox.size() || obslist.size() > numobs) ? '+' : '.' } ;
	cout << c << flush ;
	}
    inbox.clear() ;
    }

void Build::do_geo_bckt (const uint3& bckt)	// Do bucket of geodesic eqns
    {
    if (global.interrupt) return ;

    uint	blab	{ Blab::blablevel[BLAB::BUILD] } ;
    auto&	info	{ global.info() } ;
    auto&	list	{ ObsList::obs } ;
    const auto&	gens	{ info.gens.front() } ;
    ushort	ngens	{ info.neven.front() } ;
    uint	offset	{ global.stage ? global.nobsG() : 0 } ;
    uint	bcktnum	{ bckt[0] } ;
    uint	first	{ bckt[1] } ;
    uint	last	{ bckt[2] } ;
    uint	bcktsiz	{ last - first + 1 } ;
    auto&	geos	{ global.data().geos[bcktnum] } ;
    PolyMap	ans	{ ObsList::obs } ;

    ulong		ngeo	(0) ;
    ulong		nterms	(0) ;
    array<ulong,PSIZ+1>	geotermord {} ;

    geos.clear() ;
    for (uint i(first) ; i <= last ; ++i)
	{
	if (list(i).is_fermi() != global.stage)
	    fatal ("do_geo_bckt: bad obs stage") ;

	if (list(i).type == ObsType::Eloop) continue ;
	for (ushort k(0) ; k < ngens ; ++k)
	    {
	    if (gens[k].T_odd)	  continue ;
	    if (global.interrupt) break ;
	    if (blab > 4) cout << "\ndo_geo_bckt " << bcktnum << " [gen " << k
			       << ", " << list(i) << "]\n" << flush ;

	    Commute::commute_poly (gens[k], ObsPoly(i,list).negate(), ans) ;
	    if (theory.euclid && gens[k].imag && list(i).imag()) ans.negate() ;
	    if (ans.size())
		{
		check_xorder (i, gens[k], ans) ;
		for (auto& [indx,coeff] : ans)
		    {
		    if (!coeff) continue ;
		    ++geotermord[indx.order()] ;
		    ++nterms ;
		    }
		}
	    ++ngeo ;
	    geos.add (Poly {GeoHdr {i,k}}, ans) ;
	    ans.clear() ;
	    }
	}
    if (geos.entry().id != RecordID::Geos) fatal ("mk_geos_bckt: bad record ID!") ;
    geos.shrink_to_fit() ;
    geos.entry().ncol   = bcktsiz ;
    geos.entry().nrow   = ngens ;
    geos.entry().reclen = geos.size() ;

    global.count().ngeos    += ngeo ;
    global.count().geoterms += nterms ;
    global.count().geotermord += geotermord ;

    if (blab > 2)
	{
	cout << "geo bucket " << std::setw(8) << bckt[1] << ":" << std::setw(8)
	     << std::left << bckt[2] << std::right<< " done\n" << flush ;
	}
    else if (blab) cout << '.' << flush ;
    if (global.autosave) Save::write_geo_bckt (bcktnum) ;
    }

void Build::check_xorder (uint i, const Gen& g, const PolyMap& ans)
    //
    // check expectation orders
    //
    {
    const auto&	list	{ ans.obslist() } ;
    const auto&	obs	{ list(i) } ;
    short	minord	{ SHRT_MAX } ;

    for (auto& [indx,coeff] : ans)
	{
	if (coeff)
	    {
	    short xord { g.order } ;
	    for (auto& k : indx) if (k) xord += list(k).xorder ;
	    if (xord < minord) minord = xord ;
	    }
	}
    if (minord < obs.xorder)
	{
	cout << "Warning: Obs[" << i << "] = " << obs << " xord "
	     << obs.xorder << " should be " << minord << "!\n" << flush ;
	}
    }

void Build::do_geostats()			// (Re)evaluate geodesic stats
    {
    global.count().ngeos      = 0 ;
    global.count().geoterms   = 0 ;
    global.count().geotermord = 0 ;

    auto oldstage { global.stage } ;
    for (int stage(0) ; stage < 2 ; ++stage)
	{
	global.stage = stage ? Global::Fermi : Global::Gauge ;
	const auto&	geos { global.data().geos } ;
	const auto&	bckt { global.info(stage).bckt } ;

	//std::sort (bckt.begin(), bckt.end(), [geos](const uint3& a, const uint3& b)
	//    { return geos[a[0]].entry().filepos < geos[b[0]].entry().filepos ; }) ;

	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), Build::do_geostat_bckt)) ;
	}
    global.stage = oldstage ;
    }

void Build::do_geostat_bckt (const uint3& bckt)	// Evaluate geo bucket statistics
    {
    if (global.interrupt) return ;

    uint		bcktnum { bckt[0] } ;
    uint		first	{ bckt[1] } ;
    uint		last	{ bckt[2] } ;
    ulong		ngeos	(0) ;
    ulong		nterms	(0) ;
    array<ulong,PSIZ+1>	byord	{} ;
    auto&		geos	{ global.data().geos[bcktnum] } ;

    if (global.geoswap) Save::read_geo_bckt (bcktnum) ;
    for (const auto& poly : geos)
	{
	const GeoHdr& info ( poly) ;
	if (info.indx < first || info.indx > last)
	    fatal (format ("Inconsistent geo record: bucket {} indx {} first {} last {}!",
		    bcktnum, info.indx, first, last)) ; 

	if (global.interrupt) break ;
	for (const auto& term : poly)
	    {
	    ++byord [term.order()] ;
	    ++nterms ;
	    }
	++ngeos ;
	}
    if (global.geoswap) Save::read_geo_bckt (-bcktnum-1) ;
    global.count().ngeos      += ngeos ;
    global.count().geoterms   += nterms ;
    global.count().geotermord += byord ;
    }

