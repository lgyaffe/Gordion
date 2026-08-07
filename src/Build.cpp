#include "Build.h"
#include "Canon.h"
#include "Commute.h"
#include "Numerics.h"
#include "Rep.h"
#include "Save.h"
#include "Blab.h"
#include "Gripe.h"

void Build::clear_obs (int stage)		// Clear prior observables
    {
    if (stage == 0)
	{
	Canon::cache.clear () ;
	ObsList::obs.clear () ;
	global.clearpolys (0) ;
	global.clearpolys (1) ;
	}
    else if (stage && global.info(1).nobs)
	{
	cout << "Purging fermion observables" ;
	Canon::cache.purge (global.info(0).nobs) ;
	ObsList::obs.purge (global.info(0).nobs) ;
	global.clearpolys (1) ;
	}
    }

void Build::mk_obs (int target)			// Build observables
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	info	{ global.stageinfo } ;
    auto&	maxord  { global.maxord() } ;
    auto&	obslist { ObsList::obs } ;

    obslist.approx  = false ;
    ObsList::freeze = false ;
    Canon::cache.freeze = false ;
    if (global.stage == Global::Gauge)
	{
	clear_obs (1) ;
	if (target < maxord) clear_obs (0) ;
	if (info[0].nobs <= 1) obslist.obsinit (0) ;
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
	    numerics.initialize  (0,info[0].nobs,info[1].nobs) ;
	    global.clearpolys    (0) ;
	    global.close_streams (0) ;
	    global.mk_bcktlist   ()  ;
	    }
	}
    else if (theory.nf) // stage == Global::Fermi
	{
	if (target < maxord) clear_obs (1) ;
	if (!info[1].nobs) obslist.obsinit (1) ;
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
	    mk_Efermions() ;
	    }
	vector<Op>().swap (primepOps) ;
	int initfail { obslist.do_fermiinit () } ;
	if (blab)
	    {
	    cout << "Fermion initializations: " << ObsList::fermiinit.size() ;
	    if (initfail) cout << " + " << initfail << " missing loop partners" ;
	    cout << "\n" ;
	    }
	if (obslist.size() != numobs)
	    {
	    numerics.initialize  (1,info[0].nobs,info[1].nobs) ;
	    global.clearpolys    (1) ;
	    global.close_streams (1) ;
	    global.mk_bcktlist   ()  ;
	    }
	}
    Obsset().swap (newobs) ;
    Obsset().swap (ObsList::inbox) ;
    Canon::cache.freeze = true ;
    ObsList::freeze = true ;
    obslist.approx = global.approx ;
    obslist.shrink() ;

    if (blab) cout << "Total # Obs:\t" << ObsList::obs.size() << "\n" << flush ;
    if (global.info().Hterms[0].cpoly.empty()) mk_ham () ;
    }

void Build::mk_loops ()				// Build Loop's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(0).bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
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
	    Obsset().swap	(newobs) ;
	    }
	if (blab > 1 && delta_n)
	    cout << "\n" << std::left << std::setw(12) << "Loops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# gauge obs = " << global.info(0).nobs << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "Loops     (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # gauge obs = "
	     << global.info(0).nobs << "\n" << flush ;
    }

void Build::mk_Eloops ()			// Build Eloop's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(0).bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
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
	    Obsset().swap	(newobs) ;
	    }
	if (blab > 1 && delta_n)
	    cout << "\n" << std::left << std::setw(12) << "Eloops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# gauge obs = " << global.info(0).nobs << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "Eloops    (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # gauge obs = "
	     << global.info(0).nobs << "\n" << flush ;
    }

void Build::mk_EEloops ()			// Build EEloop's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(0).bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
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
	    Obsset().swap	(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "EEloops"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# gauge obs = " << global.info(0).nobs << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "EEloops   (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # gauge obs = "
	     << global.info(0).nobs << "\n" << flush ;
    }

void Build::mk_fermions ()			// Build Fermion's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(1).bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
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
	    Obsset().swap	(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "Fermions"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# fermion obs = " << global.info(1).nobs << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "Fermions  (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # fermion obs = "
	     << global.info(1).nobs << "\n" << flush ;
    }

void Build::mk_Efermions ()			// Build Efermion's
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info(1).bckt } ;
    auto&	maxord	{ global.maxord() } ;
    auto	prev	{ global.maxthread } ;
    auto&	obslist	{ ObsList::obs } ;
    auto	numobs	{ obslist.size() } ;

    for (cord = maxord/2 ; cord <= maxord - 1 ; ++cord)
	{
	if (blab > 2) cout << "mk_Efermions: maxord " << maxord
			   << " cord " << cord << "\n" << flush ;

	int	delta_n	( -obslist.size() ) ;
	global.mk_bcktlist () ;
	if (obslist.size() <= single_thread) global.maxthread = 1 ;
	else Canon::cache.freeze = true ;
	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), do_Efermion_bckt)) ;
	global.maxthread = prev ;
	delta_n += obslist.size() + newobs.size() ;

	if (newobs.size())
	    {
	    Canon::cache.load	(newobs) ;
	    obslist.insert	(newobs) ;
	    Obsset().swap	(newobs) ;
	    }
	if (blab > 1)
	    cout << "\n" << std::left << std::setw(12) << "Efermions"
		 << "(" << cord << "," << maxord-cord << "): \t"
		 << std::right << std::setw(9) << delta_n << " added, "
		 << "# fermion obs = " << global.info(1).nobs << "\n" << flush ;
	if (global.interrupt) return ;
	}
    if (blab == 1 && numobs < obslist.size())
	cout << "Efermions (" << maxord << "):  \t"
	     << std::right << std::setw(9)
	     << obslist.size() - numobs << " added, # fermion obs = "
	     << global.info(1).nobs << "\n" << flush ;
    }

void Build::mk_ham()				// Build canonical H
    {
    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    auto&	baslist { ObsList::base } ;
    auto&	obslist { ObsList::obs  } ;
    auto&	MMAobs  { global.info().MMAfile.obs } ;
    auto&	Hterms	{ global.info().Hterms } ;
    auto	nterms	( Hterms.size() ) ;
    PolyMap	ans	{ obslist } ;
    string	name	{ theory.euclid ? "Free energy" : "Hamiltonian" } ;

    if (blab) cout << name << ": " << flush ;
    for (int i(0) ; i < nterms ; ++i)
	{
	for (const auto& term : Hterms[i].poly)
	    {
	    Obs  o	{ baslist (term[0]) } ; o.canon() ;
	    numb indx	{ obslist.find (o) } ;
	    ans.add (PolyTerm (indx, term.coeff)) ;
	    if (indx) MMAobs.try_emplace (indx, o) ;
	    }
	Hterms[i].cpoly.clear() ;
	Hterms[i].cpoly.push_map (ans) ;
	ans.clear() ;
	}
    if (blab) cout << "\t\tdone\n" << flush ;
    }

void Build::mk_grad()				// Build gradient
    {
    if (global.interrupt) return ;
    if (ObsList::obs.swapped) Save::reload_obs() ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	Hterms	{ global.info().Hterms } ;
    const auto&	gens	{ global.info().gens.front() } ;
    auto	neven	{ global.info().neven.front() } ;
    auto&	grad	{ global.data().grad } ;
    auto	nterms	( Hterms.size() ) ;
    PolyMap	ans	{ ObsList::obs } ;

    if (!global.maxord()) gripe ("Make some observables first!") ;
    if (grad.entry().id != RecordID::Grad) fatal ("mk_grad: bad record ID!") ;

    grad.clear() ;
    if (blab) cout << "Gradient:    " << flush ;
    for (int i(0) ; i < nterms ; ++i)
	{
	const ObsPoly& poly { Hterms[i].cpoly } ;
	for (int j(0) ; j < neven ; ++j)
	    {
	    Commute::commute_poly (gens[j], poly, ans) ;
	    if (Hterms[i].imag && gens[j].imag) ans.negate() ;
	    grad.add (ans.negate()) ;
	    ans.clear() ;
	    if (global.interrupt) return ;
	    }
	}
    grad.shrink_to_fit() ;
    grad.entry().ncol   = nterms ;
    grad.entry().nrow   = neven ;
    grad.entry().reclen = grad.size() ;
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
    if (ObsList::obs.swapped) Save::reload_obs() ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	repnam	{ Rep::list[repnum].name } ;
    const auto&	gens	{ global.info().gens[repnum] } ;
    const auto&	Hterms	{ global.info().Hterms } ;
    auto&	curv	{ global.data().curv[repnum] } ;
    auto	nterms	( Hterms.size() ) ;
    auto	ngens   ( gens.size() ) ;
    ObsList	tmplist { "CurvTemp", repnum == 0 } ;
    PolyMap	ans	{ ObsList::obs } ;
    PolyMap	tmp	{ tmplist } ;

    if (!global.maxord()) gripe ("Make some observables first!") ;
    if (curv.entry().id != RecordID::Curv) fatal ("mk_curv: bad record ID!") ;

    curv.clear() ;
    if (blab && ngens && nterms) cout << repnam << " curvature:    " << flush ;
    for (int i(0) ; i < nterms ; ++i)
	{
	const ObsPoly& poly { repnum ? Hterms[i].poly : Hterms[i].cpoly } ;

	for (int j(0) ; j < ngens ; ++j)
	    {
	    for (int k(0) ; k < ngens ; ++k)
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
		curv.add (ans) ;
		ans.clear() ;
		tmp.clear() ;
		if (global.interrupt) return ;
		}
	    }
	}
    curv.shrink_to_fit() ;
    curv.entry().nslice  = nterms ;
    curv.entry().ncol    = ngens ;
    curv.entry().nrow    = ngens ;
    curv.entry().reclen  = curv.size() ;
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
    if (ObsList::obs.swapped) Save::reload_obs() ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	repnam	{ Rep::list[repnum].name } ;
    const auto&	gens	{ global.info().gens[repnum] } ;
    auto	neven	{ global.info().neven[repnum] } ;
    auto&	lagr	{ global.data().lagr[repnum] } ;
    auto&	oplist	{ global.info().ops } ;
    auto	opnum	{ oplist.size() } ;
    auto	ngens	( gens.size() ) ;
    short	trunc	{ SHRT_MAX } ;
    PolyMap	ans	{ ObsList::obs } ;

    if (!global.maxord()) gripe ("Make some observables first!") ;
    if (lagr.entry().id != RecordID::Lagr) fatal ("mk_lagr: bad record ID!") ;

    lagr.clear();
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
    lagr.shrink_to_fit() ;
    lagr.entry().ncol   = neven ;
    lagr.entry().nrow   = ngens - neven ;
    lagr.entry().reclen = lagr.size();
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
    if (ObsList::obs.swapped) Save::reload_obs() ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	bckt	{ global.info().bckt } ;

    if (blab && blab < 3) cout << "geodesics: " << flush;
    if (!global.maxord()) gripe ("Make some observables first!") ;

    ObsList::freeze = true ;
    Canon::cache.freeze = true ;
    global.count().cleargeostats() ;
    if (global.autosave) Save::save_sys() ;
    TASK_ARENA (global.maxthread, bckt,
	FOR_EACH (bckt.begin(), bckt.end(), do_geo_bckt)) ;

    if (blab) cout << "\tdone\n" << flush ;
    }

void Build::do_Loop_bckt (const numb3& bckt)		// Do loop build bucket 
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	oplist	{ global.info(0).ops } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    auto&	obslist	{ ObsList::obs } ;
    PolyMap 	tmp	{ obslist } ;
    PolyTerm 	zero	{ PolyIndx(), 0 } ;
    auto	numobs	{ obslist.size() } ;
    numb	bcktnum { bckt[0] } ;
    numb 	first	{ bckt[1] } ;
    numb 	last	{ bckt[2] } ;

    ObsList::freeze = false ;
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

void Build::do_Eloop_bckt (const numb3& bckt)		// Do Eloop build bucket
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	oplist	{ global.info(0).ops } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    auto&	obslist	{ ObsList::obs } ;
    PolyTerm	zero	{ PolyIndx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    auto	numobs	{ obslist.size() } ;
    numb	bcktnum { bckt[0] } ;
    numb 	first	{ bckt[1] } ;
    numb 	last	{ bckt[2] } ;

    ObsList::freeze = false ;
    for (numb i(first) ; i <= last ; ++i)
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

void Build::do_EEloop_bckt (const numb3& bckt)		// Do EEloop build bckt
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	oplist	{ global.info(0).ops } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    auto&	obslist	{ ObsList::obs } ;
    auto	numobs	{ obslist.size() } ;
    PolyTerm	zero	{ PolyIndx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    numb	bcktnum { bckt[0] } ;
    numb 	first	{ bckt[1] } ;
    numb 	last	{ bckt[2] } ;

    ObsList::freeze = false ;
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

void Build::do_Fermion_bckt (const numb3& bckt)		// Do Fermion build bckt
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    auto&	obslist	{ ObsList::obs } ;
    auto	numobs	{ obslist.size() } ;
    PolyTerm	zero	{ PolyIndx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    numb	bcktnum { bckt[0] } ;
    numb	first	{ bckt[1] } ;
    numb	last	{ bckt[2] } ;

    ObsList::freeze = false ;
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

void Build::do_Efermion_bckt (const numb3& bckt)	// Do Efermion build bckt
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	oplist	{ global.info(1).ops } ;
    const auto&	maxord	{ global.maxord() } ;
    auto&	inbox	{ ObsList::inbox } ;
    auto&	obslist	{ ObsList::obs } ;
    auto	numobs	{ obslist.size() } ;
    PolyTerm	zero	{ PolyIndx(), 0 } ;
    PolyMap	tmp	{ obslist } ;
    numb	bcktnum { bckt[0] } ;
    numb 	first	{ bckt[1] } ;
    numb 	last	{ bckt[2] } ;

    ObsList::freeze = false ;
    for (numb i(first) ; i <= last ; ++i)
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
		if (global.interrupt) return ;
		}
	    }
	else if (obslist(i).is_Efermion())
	    {
	    Obs a	{ obslist(i) } ;
	    int ord	{ a.order() } ;
	    for (const auto& op : primepOps)
		{
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
		if (global.interrupt) return ;
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

void Build::do_geo_bckt (const numb3& bckt)		// Do bucket of geodesic eqns
    {
    if (global.interrupt) return ;

    const auto&	blab	{ Blab::level(Blab::BUILD) } ;
    const auto&	gens	{ global.info().gens.front() } ;
    auto	neven	{ global.info().neven.front() } ;
    auto&	list	{ ObsList::obs } ;
    numb	bcktnum	{ bckt[0] } ;
    numb	first	{ bckt[1] } ;
    numb	last	{ bckt[2] } ;
    numb	bcktsiz	{ last - first + 1 } ;
    auto&	geos	{ global.data().geos[bcktnum] } ;
    PolyMap	ans	{ ObsList::obs } ;

    ulong		ngeo	(0) ;
    ulong		nterms	(0) ;
    array<ulong,PSIZ+1>	geotermord {} ;

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

    global.count().ngeos      += ngeo ;
    global.count().geoterms   += nterms ;
    global.count().geotermord += geotermord ;

    if (blab > 2)
	{
	cout << "geo bucket " << std::setw(8) << bckt[1] << ":" << std::setw(8)
	     << std::left << bckt[2] << std::right<< " done\n" << flush ;
	}
    else if (blab) cout << '.' << flush ;
    if (global.autosave) Save::write_geo_bckt (bcktnum) ;
    }

void Build::check_xorder (numb i, const Gen& g, const PolyMap& ans)	// Check xorders
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
    for (int stage(0) ; stage < 2 - !theory.nf ; ++stage)
	{
	global.stage = stage ? Global::Fermi : Global::Gauge ;
	const auto&	geos { global.data().geos } ;
	const auto&	bckt { global.info(stage).bckt } ;

	TASK_ARENA (global.maxthread, bckt,
	    FOR_EACH (bckt.begin(), bckt.end(), Build::do_geostat_bckt)) ;
	}
    global.stage = oldstage ;
    }

void Build::do_geostat_bckt (const numb3& bckt)	// Evaluate geo bucket statistics
    {
    if (global.interrupt) return ;

    const auto&		stage	{ global.stage } ;
    const auto&		bcktnum { bckt[0] } ;
    const auto&		first	{ bckt[1] } ;
    const auto&		last	{ bckt[2] } ;
    const auto&		geos	{ global.data().geos[bcktnum] } ;
    ulong		ngeos	(0) ;
    ulong		nterms	(0) ;
    array<ulong,PSIZ+1>	byord	{} ;

    if (global.geoswap) Save::read_geo_bckt (stage, bcktnum) ;
    const auto&		ncol { geos.entry().ncol } ;
    const auto&		nrow { geos.entry().nrow } ;

    for (const auto& poly : geos)
	{
	if (global.interrupt) return ;
	auto ptr { poly.begin() } ;
	for (const auto& term : poly)
	    {
	    PolyTerm t { poly.nextterm (ptr) } ;
	    ++byord [t.order()] ;
	    ++nterms ;
	    }
	++ngeos ;

	if (geos.entry().items() != ngeos)
	    fatal (format("Inconsistent geo record: bucket {} expected {} got {}",
		    bcktnum, geos.entry().items(), ngeos)) ;
	}
    if (global.geoswap) Save::read_geo_bckt (stage, -bcktnum-1) ;
    global.count().ngeos      += ngeos ;
    global.count().geoterms   += nterms ;
    global.count().geotermord += byord ;
    }

