#include "Build.h"
#include "Canon.h"
#include "Commute.h"
#include "Global.h"
#include "Parallel.h"
#include "Rep.h"
#include "Save.h"
#include "Blab.h"
#include "Gripe.h"
#include <mutex>

static uint		cord ;			// Current creation order
static ObsSet	 	newobs ;		// Newly generated Obs
static vector<Op>	primepOps ;		// primary conj mom Ops
static std::mutex	obsmutex ;		// newobs insertion lock
static constexpr int	single_thread { 1024 };	// multithreading threshold

void Build::mk_eqns (string word)		// Build polynomial scripts
    {
    try { mk_eqns (Rep::known (word)) ; }
    catch (const exception& e) { gripe ("Unknown representation " + word) ; }
    }

void Build::mk_eqns (uint rep)			// Build polynomial scripts
    {
    bool	isH { !theory.euclid } ;
    if (!rep) { mk_ham  () ;
		mk_grad () ; }
		mk_curv (rep) ;
    if (isH)	mk_lagr (rep) ;
    if (!rep)	mk_geos () ;
    }

void Build::mk_ham()				// Build canonical H
    {
    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	maxord	{ global.info().maxord } ;
    ObsList&	baslist { global.baseobs } ;
    ObsList&	obslist { global.obs  } ;
    auto&	MMAobs  { global.info().MMAfile.obs } ;
    auto&	Hterms	{ global.info().Hterms } ;
    auto&	ham	{ global.data().ham } ;
    auto	nterms	{ Hterms.size() } ;
    PolyMap	ans	{ obslist } ;
    string	label	{ theory.euclid ? "Free energy" : "Hamiltonian" } ;

    if (ham.entry().id != RecordID::Ham) fatal ("mk_ham: bad record ID!") ;
    if (!maxord && !global.xtraobs) gripe ("Make some observables first!") ;

    global.obs.freeze = false ;
    ham.clear() ;
    if (blab) cout << label << ": " << flush ;
    for (int i(0) ; i < nterms ; ++i)
	{
	for (const auto& Hterm : Hterms[i].poly)
	    {
	    Obs  obs	{ baslist (Hterm[0]) } ; obs.canon() ;
	    auto term	{ obslist.catalog (obs) } ;
	    auto indx	{ term.item[0] } ;
	    ans.add (term * Hterm.coeff) ;
	    if (indx) MMAobs.try_emplace (obs, indx) ;
	    }
	ham.add (ans) ;
	ans.clear() ;
	}
    ham.shrink_to_fit () ;
    ham.entry().nelem  = nterms ;
    ham.entry().reclen = ham.size() ;
    global.obs.freeze  = true ;
    if (blab) cout << "\t\tdone\n" << flush ;
    }

void Build::mk_grad()				// Build gradient
    {
    if (global.interrupt) return ;
    if (global.obs.swapped) Save::reload_obs() ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	gens	{ global.info().gens.front() } ;
    const auto&	Hterms	{ global.info().Hterms } ;
    const auto& maxord	{ global.info().maxord } ;
    const auto&	ham	{ global.data().ham } ;
    auto&	grad	{ global.data().grad } ;
    auto& 	stats	{ global.info().count.gradstats } ;
    auto	nterms	{ Hterms.size() } ;

    if (ham.entry().nelem != nterms) mk_ham() ;
    if (grad.entry().id != RecordID::Grad) fatal ("mk_grad: bad record ID!") ;
    if (!maxord && !global.xtraobs) gripe ("Make some observables first!") ;
    if (!gens.size()) Gen::geninit (global.stage) ;

    const auto&	neven	{ global.info().neven.front() } ;
    const auto&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;
    PolyMap	ans	{ global.obs } ;

    global.obs.freeze = !global.xtraobs ;
    grad.clear() ;
    stats.clear() ;
    if (blab) cout << "Gradient:    " << flush ;
    for (int i(0) ; i < nterms ; ++i)
	{
	const ObsPoly poly { global.obs, ham(i) } ;
	for (int j(0) ; j < neven ; ++j)
	    {
	    Commute::commute_poly (gens[j], poly, ans) ;
	    if (Hterms[i].imag && gens[j].imag) ans.negate() ;
	    grad.add (ans.negate()) ;
	    ans.clear() ;
	    if (global.interrupt) return ;
	    }
	}
    if (obslist.size() != numobs) added_obs () ;
    grad.shrink_to_fit() ;
    grad.entry().ncol   = nterms ;
    grad.entry().nrow   = neven ;
    grad.entry().reclen = grad.size() ;
    global.obs.freeze	= true ;
    stats.add (grad) ;
    if (blab) cout << "\t\tdone\n" << flush ;
    }

void Build::mk_curv (string word)			// Build curvature
    {
    try { mk_curv (Rep::known (word)) ; }
    catch (const exception& e) { gripe ("Unknown representation " + word) ; }
    }

void Build::mk_curv (uint repnum)			// Build curvature
    {
    if (global.interrupt) return ;
    if (global.obs.swapped) Save::reload_obs() ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	repnam	{ Rep::list[repnum].name } ;
    const auto&	obslist	{ global.obs } ;
    const auto&	gens	{ global.info().gens[repnum] } ;
    const auto&	Hterms	{ global.info().Hterms } ;
    const auto& maxord	{ global.info().maxord } ;
    const auto&	ham	{ global.data().ham } ;
    auto&	curv	{ global.data().curv[repnum] } ;
    auto& 	stats	{ global.info().count.curvstats } ;
    auto	nterms	{ Hterms.size() } ;

    if (ham.entry().nelem != nterms) mk_ham() ;
    if (curv.entry().id != RecordID::Curv) fatal ("mk_curv: bad record ID!") ;
    if (!maxord && !global.xtraobs) gripe ("Make some observables first!") ;
    if (!global.info().gens[0].size()) Gen::geninit (global.stage) ;

    auto	numobs	{ obslist.size() } ;
    auto	ngens	{ gens.size() } ;
    ObsList	tmplist	{ "CurvTemp", repnum == 0 } ;
    PolyMap	ans	{ global.obs } ;
    PolyMap	tmp	{ tmplist } ;

    global.obs.freeze = !global.xtraobs ;
    curv.clear() ;
    stats.clear() ;
    if (blab && ngens && nterms) cout << repnam << " curvature:    " << flush ;
    for (int i(0) ; i < nterms ; ++i)
	{
	const ObsPoly	poly { global.obs, ham(i) } ;
	const ObsPoly*	ptr  { repnum ? &Hterms[i].poly : &poly } ;

	for (int j(0) ; j < ngens ; ++j)
	    {
	    for (int k(0) ; k < ngens ; ++k)
		{
		if (gens[j].T_odd == gens[k].T_odd)
		    {
		    Commute::commute_poly (gens[j], *ptr, tmp) ;
		    Commute::commute_poly (gens[k], tmp,  ans) ;
		    if (gens[j].imag != gens[k].imag && !theory.euclid)
			fatal ("mk_curv: imaginary!") ;
		    else if (gens[j].imag && gens[k].imag   ||
			     gens[j].imag && Hterms[i].imag ||
			     gens[k].imag && Hterms[i].imag) ans.negate() ;
		    }
		curv.add (ans) ;
		ans.clear() ;
		tmp.clear() ;
		if (global.interrupt) return ;
		}
	    }
	}
    if (obslist.size() != numobs) added_obs () ;
    curv.shrink_to_fit () ;
    curv.entry().nslice = nterms ;
    curv.entry().ncol   = ngens ;
    curv.entry().nrow   = ngens ;
    curv.entry().reclen = curv.size() ;
    global.obs.freeze	= true ;
    stats.add (curv) ;
    if (blab && ngens && nterms) cout << "\tdone\n" << flush ;
    }

void Build::mk_lagr (string word)			// Build Lagrange bracket
    {
    try { mk_lagr (Rep::known (word)) ; }
    catch (const exception& e) { gripe ("Unknown representation " + word) ; }
    }

void Build::mk_lagr (uint repnum)			// Build Lagrange bracket
    {
    if (global.interrupt) return ;
    if (global.obs.swapped) Save::reload_obs() ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	repnam	{ Rep::list[repnum].name } ;
    const auto&	obslist	{ global.obs } ;
    const auto& maxord	{ global.info().maxord } ;
    const auto&	gens	{ global.info().gens[repnum] } ;
    auto&	lagr	{ global.data().lagr[repnum] } ;
    auto& 	stats	{ global.info().count.lagrstats } ;

    if (!maxord && !global.xtraobs) gripe ("Make some observables first!") ;
    if (lagr.entry().id != RecordID::Lagr) fatal ("mk_lagr: bad record ID!") ;
    if (!global.info().gens[0].size()) Gen::geninit (global.stage) ;

    const auto&	neven	{ global.info().neven[repnum] } ;
    auto	numobs	{ obslist.size() } ;
    auto&	oplist	{ global.info().ops } ;
    auto	opnum	{ oplist.size() } ;
    short	trunc	{ SHRT_MAX } ;
    PolyMap	ans	{ global.obs } ;

    auto ngens { gens.size() } ;
    global.obs.freeze = !global.xtraobs ;
    lagr.clear();
    stats.clear() ;
    if (ngens && blab) cout << repnam << " Lagrange brkt: " << flush ;
    for (int j(0) ; j < neven ; ++j)
	{
	for (int k(neven) ; k < ngens ; ++k)
	    {
	    Gen newgen { oplist } ;
	    Commute::commute_gen (gens[j], gens[k], newgen) ;
	    if (!ans.add_gen (newgen))
		trunc = std::min(trunc, newgen.order) ;
	    lagr.add (ans) ;
	    ans.clear() ;
	    if (global.interrupt) return ;
	    }
	}
    if (obslist.size() != numobs) added_obs () ;
    lagr.shrink_to_fit() ;
    lagr.entry().ncol   = neven ;
    lagr.entry().nrow   = ngens - neven ;
    lagr.entry().reclen = lagr.size();
    global.obs.freeze	= true ;
    stats.add (lagr) ;
    if (blab && ngens)
	{
	if (trunc == SHRT_MAX) cout << "\tdone\n" ;
	else cout << format("\ttruncation in creation order {} elements\n", trunc) ;
	}
    oplist.purge (opnum) ;
    }

void Build::mk_geos()				// Build geodesic equations
    {
    if (global.interrupt) return ;
    if (global.obs.swapped) Save::reload_obs() ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info().bckt } ;
    const auto&	gens	{ global.info().gens.front() } ;
    auto&	stats	{ global.info().count.geostats } ;

    if (blab && blab < 3) cout << "geodesics: " << flush;
    if (!global.info().maxord) gripe ("Make some observables first!") ;
    if (!gens.size()) Gen::geninit (global.stage) ;

    global.obs.freeze   = true ;
    Canon::cache.freeze = true ;
    stats.clear() ;
    if (Save::autosave) Save::save_sys() ;
    TASK_ARENA (global.maxthread, bckt,
	FOR_EACH (bckt.begin(), bckt.end(), do_geo_bckt)) ;

    if (blab) cout << "\tdone\n" << flush ;
    }

void Build::mk_obs (uint target)		// Build observables
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	info	{ global.stageinfo } ;
    const auto&	gens	{ global.info().gens.front() } ;
    auto&	maxord  { global.info().maxord } ;
    ObsList&	obslist { global.obs } ;

    obslist.classify	= true ;
    obslist.approx	= false ;
    global.obs.freeze	= false ;
    Canon::cache.freeze = false ;
    if (global.stage == Global::Gauge)
	{
	clear_obs (1) ;
	if (target < maxord) clear_obs (0) ;
	if (!gens.size()) Gen::geninit (0) ;
	if (obslist.nobsG <= 1) obslist.obsinit (0) ;
	auto numobs { obslist.size() } ;
	while (maxord < target)
	    {
	    if (global.interrupt) return ;
	    if (++maxord % 4) continue ;
	    mk_loops () ;
	    if (theory.euclid) continue ;
	    mk_Eloops() ;
	    mk_EEloops() ;
	    }
	if (obslist.size() != numobs)
	    {
	    added_obs () ;
	    global.clearpolys (0) ;
	    }
	}
    else if (theory.nf) // stage == Global::Fermi
	{
	if (target < maxord) clear_obs (1) ;
	if (!obslist.nobsF) obslist.obsinit (1) ;
	auto& oplistG { info[0].ops } ;
	auto& oplistF { info[1].ops } ;
	auto  numobs  { obslist.size() } ;
	auto  lambda  { [](const Op& o) { return o.primary && !o.is_coord(); } } ;

	std::copy_if (oplistG.begin(), oplistG.end(),
		      std::back_inserter (primepOps), lambda) ;
	std::copy_if (oplistF.begin(), oplistF.end(),
		      std::back_inserter (primepOps), lambda) ;

	while (maxord < target)
	    {
	    if (global.interrupt) return ;
	    if (++maxord % 2) continue ;
	    mk_fermions () ;
	    if (theory.euclid) continue ;
	    mk_efermions() ;
	    }
	vector<Op>().swap (primepOps) ;
	int initfail { obslist.do_fermiinit () } ;
	if (blab)
	    {
	    cout << "Fermion initializations: " << global.obs.fermiinit.size() ;
	    if (initfail) cout << " + " << initfail << " missing loop partners" ;
	    cout << "\n" ;
	    }
	if (obslist.size() != numobs)
	    {
	    added_obs () ;
	    global.clearpolys (1) ;
	    }
	}
    ObsSet().swap (newobs) ;
    ObsSet().swap (global.obs.inbox) ;
    Canon::cache.freeze	= true ;
    global.obs.freeze	= true ;
    obslist.classify	= false ;
    obslist.approx	= global.approx ;
    obslist.shrink() ;

    if (blab) cout << "Total # Obs:\t" << global.obs.size() << "\n" << flush ;
    }

void Build::mk_loops ()				// Build Loop's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(0).bckt } ;
    auto&	maxord	{ global.info(0).maxord } ;
    auto	prev	{ global.maxthread } ;
    ObsList&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;

    for (cord = maxord/2 ; cord <= maxord - 4 ; cord += 2)
	{
	if (blab > 2) cout << "\nmk_loops: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.size() ) ;
	global.mk_bcktlist () ;
	if (obslist.size() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_Loop_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.size() + newobs.size() ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    ObsSet().swap	(newobs) ;
	    }
	if (blab > 1 && delta_n)
	    cout << "\n" << std::left << std::setw(12) << "Loops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# gauge obs = " << obslist.nobsG << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "Loops     (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # gauge obs = "
	     << obslist.nobsG << "\n" << flush ;
    }

void Build::mk_Eloops ()			// Build Eloop's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(0).bckt } ;
    auto&	maxord	{ global.info(0).maxord } ;
    auto	prev	{ global.maxthread } ;
    ObsList&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;
    int		maxcord	{ maxord - (maxord > 4 ? 4 : 2) } ;

    for (cord = maxord/2 ; cord <= maxcord ; cord += 2)
	{
	if (blab > 2) cout << "\nmk_Eloops: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.size() ) ;
	global.mk_bcktlist () ;
	if (obslist.size() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_Eloop_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.size() + newobs.size() ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    ObsSet().swap	(newobs) ;
	    }
	if (blab > 1 && delta_n)
	    cout << "\n" << std::left << std::setw(12) << "Eloops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# gauge obs = " << obslist.nobsG << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "Eloops    (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # gauge obs = "
	     << obslist.nobsG << "\n" << flush ;
    }

void Build::mk_EEloops ()			// Build EEloop's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(0).bckt } ;
    auto&	maxord	{ global.info(0).maxord } ;
    auto	prev	{ global.maxthread } ;
    ObsList&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;

    for (cord = maxord/2 - 2 ; cord <= maxord - 4 ; cord += 2)
	{
	if (blab > 2) cout << "\nmk_EEloops: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.size() ) ;
	global.mk_bcktlist () ;
	if (obslist.size() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_EEloop_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.size() + newobs.size() ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    ObsSet().swap	(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "EEloops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# gauge obs = " << obslist.nobsG << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "EEloops   (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # gauge obs = "
	     << obslist.nobsG << "\n" << flush ;
    }

void Build::mk_fermions ()			// Build Fermion's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(1).bckt } ;
    auto&	maxord	{ global.info(1).maxord } ;
    auto	prev	{ global.maxthread } ;
    ObsList&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;

    for (cord = maxord/2 ; cord <= maxord - 2 ; ++cord)
	{
	if (blab > 2) cout << "mk_fermions: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.size() ) ;
	global.mk_bcktlist () ;
	if (obslist.size() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_Fermion_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.size() + newobs.size() ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    ObsSet().swap	(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "Fermions"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# fermion obs = " << obslist.nobsF << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "Fermions  (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # fermion obs = "
	     << obslist.nobsF << "\n" << flush ;
    }

void Build::mk_efermions ()			// Build efermion's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(1).bckt } ;
    auto&	maxord	{ global.info(1).maxord } ;
    auto	prev	{ global.maxthread } ;
    ObsList&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;

    for (cord = maxord/2 ; cord <= maxord - 1 ; ++cord)
	{
	if (blab > 2) cout << "mk_efermions: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.size() ) ;
	global.mk_bcktlist () ;
	if (obslist.size() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_efermion_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.size() + newobs.size() ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    ObsSet().swap	(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "efermions"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# fermion obs = " << obslist.nobsF << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "efermions (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # fermion obs = "
	     << obslist.nobsF << "\n" << flush ;
    }

void Build::do_Loop_bckt (const numb3& bckt)		// Do loop build bucket 
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	oplist	{ global.info(0).ops } ;
    const auto&	maxord	{ global.info(0).maxord } ;
    ObsList&	obslist	{ global.obs } ;
    auto&	inbox	{ global.obs.inbox } ;
    PolyMap 	tmp	{ obslist } ;
    PolyTerm 	zero	{ PolyIndx(), 0 } ;
    auto	numobs	{ obslist.size() } ;
    numb	bcktnum { bckt[0] } ;
    numb 	first	{ bckt[1] } ;
    numb 	last	{ bckt[2] } ;

    global.obs.freeze = false ;
    for (numb i(first) ; i <= last ; ++i)
	{
	if (obslist(i).type != ObsType::Loop)	continue ;

	Obs a	{ obslist(i) } ;
	int ord	{ a.order() } ;
	for (const auto& op : oplist)
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
	    if (global.interrupt) break ;
	    }
	}
    global.obs.freeze = true ;
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

void Build::do_Eloop_bckt (const numb3& bckt)		// Do Eloop build bucket
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	oplist	{ global.info(0).ops } ;
    const auto&	maxord	{ global.info(0).maxord } ;
    auto&	inbox	{ global.obs.inbox } ;
    ObsList&	obslist	{ global.obs } ;
    PolyTerm	zero	{ PolyIndx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    auto	numobs	{ obslist.size() } ;
    numb	bcktnum { bckt[0] } ;
    numb 	first	{ bckt[1] } ;
    numb 	last	{ bckt[2] } ;

    global.obs.freeze = false ;
    for (numb i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	auto	insiz { inbox.size() } ;
	if (obslist(i).is_Loop())
	    {
	    if (obslist(i).corder != cord) continue ;

	    Op op { obslist(i) } ;
	    for (const Obs* EE : global.baseobs)
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
		if (global.interrupt) return ;
		}
	    }
	else if (obslist(i).is_Eloop())
	    {
	    Obs	a	{ obslist(i) } ;
	    int	ord	{ a.order() } ;
	    for (const auto& op : oplist)
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
		if (global.interrupt) return ;
		}
	    }
	}
    global.obs.freeze = true ;
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

void Build::do_EEloop_bckt (const numb3& bckt)		// Do EEloop build bckt
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	oplist	{ global.info(0).ops } ;
    const auto&	maxord	{ global.info(0).maxord } ;
    auto&	inbox	{ global.obs.inbox } ;
    ObsList&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;
    PolyTerm	zero	{ PolyIndx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    numb	bcktnum { bckt[0] } ;
    numb 	first	{ bckt[1] } ;
    numb 	last	{ bckt[2] } ;

    global.obs.freeze = false ;
    for (numb i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	if (obslist(i).is_fermi() || !obslist(i).is_EEloop()) continue ;

	Obs a	{ obslist(i) } ;
	int ord	{ a.order() } ;
	for (const auto& op : oplist)
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
	    if (global.interrupt) return ;
	    }
	}
    global.obs.freeze = true ;
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

void Build::do_Fermion_bckt (const numb3& bckt)		// Do Fermion build bckt
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	maxord	{ global.info(1).maxord } ;
    auto&	inbox	{ global.obs.inbox } ;
    ObsList&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;
    PolyTerm	zero	{ PolyIndx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    numb	bcktnum { bckt[0] } ;
    numb	first	{ bckt[1] } ;
    numb	last	{ bckt[2] } ;

    global.obs.freeze = false ;
    for (numb i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	if (obslist(i).type != ObsType::Fermion) continue ;

	Obs a { obslist(i) } ;
	int ord	{ a.order() } ;
	for (const auto& op : primepOps)
	    {
	    if (op.order + a.corder != cord)	continue ;
	    if (ord < maxord - 2 * op.order)	continue ;
	    if (blab > 3)
		{
		cout << "do_Fermion_bckt " << bcktnum
		     << " maxord " << maxord
		     << " cord " << cord 
		     << " [" << op << "," << a
		     << "] (" << op.order << "+" << a.corder << ")\n" << flush ;
		}
	    Commute::do_commute (op, a, zero, obslist, tmp) ;
	    if (global.interrupt) return ;
	    }
	}
    global.obs.freeze = true ;
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

void Build::do_efermion_bckt (const numb3& bckt)	// Do efermion build bckt
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	oplist	{ global.info(1).ops } ;
    const auto&	maxord	{ global.info(1).maxord } ;
    auto&	inbox	{ global.obs.inbox } ;
    ObsList&	obslist	{ global.obs } ;
    auto	numobs	{ obslist.size() } ;
    PolyTerm	zero	{ PolyIndx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    numb	bcktnum { bckt[0] } ;
    numb 	first	{ bckt[1] } ;
    numb 	last	{ bckt[2] } ;

    global.obs.freeze = false ;
    for (numb i(first) ; i <= last ; ++i)
	{
	if (global.interrupt) break ;
	if (obslist(i).is_Fermion())
	    {
	    if (obslist(i).corder != cord) continue ;

	    Op op { obslist(i) } ;
	    if (!theory.euclid) op.front() = stag(op.front()) ;
	    for (const Obs* EE : global.baseobs)
		{
		if (EE->size() != 1 || !isEE(EE->front())) continue ;
		if (blab > 3)
		    {
		    cout << "\nmk_efermions " << bcktnum << " [" << op << "," << *EE
			 << "] (" << op.order << "+" << EE->corder << ")\n" << flush ;
		    }
		Commute::do_commute (op, *EE, zero, obslist, tmp) ;
		if (global.interrupt) return ;
		}
	    }
	else if (obslist(i).is_efermion())
	    {
	    Obs a	{ obslist(i) } ;
	    int ord	{ a.order() } ;
	    for (const auto& op : primepOps)
		{
		if (op.order + a.corder != cord)	continue ;
		if (ord < maxord - 2 * op.order)	continue ;
		if (blab > 3)
		    {
		    cout << "\nmk_efermions [" << op << "," << a << "] ("
			 << op.order << "+" << a.corder << ")\n" << flush ;
		    }
		if (op.is_Eloop())
		    Commute::do_commute  (op, a, zero, obslist, tmp) ;
		else
		    Commute::do_commuteE (op, a, zero, obslist, tmp) ;
		if (global.interrupt) return ;
		}
	    }
	}
    global.obs.freeze = true ;
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

void Build::do_geo_bckt (const numb3& bckt)		// Do bucket of geodesic eqns
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	gens	{ global.info().gens.front() } ;
    const auto&	neven	{ global.info().neven.front() } ;
    ObsList&	list	{ global.obs } ;
    numb	bcktnum	{ bckt[0] } ;
    numb	first	{ bckt[1] } ;
    numb	last	{ bckt[2] } ;
    numb	bcktsiz	{ last - first + 1 } ;
    auto&	geos	{ global.data().geos[bcktnum] } ;
    auto&	stats	{ global.info().count.geostats } ;
    PolyMap	ans	{ global.obs } ;

    if (geos.entry().id != RecordID::Geos) fatal ("do_geos_bckt: bad record ID!") ;
    geos.clear() ;
    for (numb i(first) ; i <= last ; ++i)
	{
	if (list(i).is_fermi() != global.stage)
	    fatal ("do_geo_bckt: bad obs stage") ;

	for (int k(0) ; k < neven ; ++k)
	    {
	    if (global.interrupt) break ;
	    if (list(i).type != ObsType::Eloop)
		{
		if (blab > 4) cout << "\ndo_geo_bckt " << bcktnum << " [gen " << k
				   << ", " << list(i) << "]\n" << flush ;

		Commute::commute_poly (gens[k], ObsPoly(list,i).negate(), ans) ;
		if (theory.euclid && gens[k].imag && list(i).imag()) ans.negate() ;
		if (ans.size())
		    {
		    check_xorder (i, gens[k], ans) ;
		    }
		}
	    geos.add (ans) ;
	    ans.clear() ;
	    if (global.interrupt) return ;
	    }
	}
    geos.shrink_to_fit() ;
    geos.entry().ncol   = bcktsiz ;
    geos.entry().nrow   = neven ;
    geos.entry().reclen = geos.size() ;
    stats.add (geos) ;

    if (blab > 2)
	{
	cout << "geo bucket " << std::setw(8) << bckt[1] << ":" << std::setw(8)
	     << std::left << bckt[2] << std::right<< " done\n" << flush ;
	}
    else if (blab) cout << '.' << flush ;
    if (Save::autosave) Save::write_geo_bckt (bcktnum) ;
    }

void Build::do_geostats()			// Reevaluate geodesic stats
    {
    auto oldstage { global.stage } ;

    for (uint stage(0) ; stage < 2 - !theory.nf ; ++stage)
	{
	global.stage = stage ? Global::Fermi : Global::Gauge ;
	const auto&	geos  { global.data().geos } ;
	const auto&	bckt  { global.info(stage).bckt } ;
	auto&		stats { global.info(stage).count.geostats } ;

	stats.clear() ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), Build::do_geostat_bckt)) ;
	}
    global.stage = oldstage ;
    }

void Build::do_geostat_bckt (const numb3& bckt) // Reevaluate geo bucket stats
    {
    if (global.interrupt) return ;

    const auto&	bcktnum { bckt[0] } ;
    const auto&	stage	{ global.stage } ;
    const auto&	geos	{ global.data().geos[bcktnum] } ;
    auto&	stats	{ global.info().count.geostats } ;

    if (global.geoswap) Save::read_geo_bckt (stage, bcktnum) ;
    stats.add (geos) ;
    if (global.geoswap) Save::read_geo_bckt (stage, -bcktnum-1) ;
    }

void Build::check_xorder (numb i, const Gen& g, const PolyMap& ans) // Check xorders
    {
    if (global.interrupt) return ;

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
	if (global.xtraobs) obs.xorder = minord ;
	else cout << "Warning: Obs[" << i << "] = " << obs << " xord "
		  << obs.xorder << " should be " << minord << "!\n" << flush ;
	}
    }

void Build::clear_obs (uint stage)		// Clear prior observables
    {
    if (stage == 0)
	{
	Canon::cache.clear () ;
	global.obs.clear  () ;
	global.clearpolys (0) ;
	global.clearpolys (1) ;
	}
    else if (stage && global.obs.nobsF)
	{
	cout << "Purging fermion observables" ;
	Canon::cache.purge (global.obs.nobsG) ;
	global.obs.purge   (global.obs.nobsG) ;
	global.clearpolys (1) ;
	}
    }

void Build::added_obs ()			// Update obs hashes
    {
    ObsList& obslist { global.obs } ;
    if (global.stage == Global::Gauge)
	{
	global.info(0).validvev = false ;
	global.info(1).validvev = false ;
	global.obs.nobsG = obslist.nobsG ;
	global.obs.nobsF = obslist.nobsF ;
	global.obs.hashG = obslist.hashG ;
	global.obs.hashF = obslist.hashF ;
	global.close_streams (0) ;
	global.mk_bcktlist   ()  ;
	}
    else
	{
	global.info(1).validvev = false ;
	global.obs.nobsF = obslist.nobsF ;
	global.obs.hashF = obslist.hashF ;
	global.close_streams (1) ;
	global.mk_bcktlist   ()  ;
	}
    }

