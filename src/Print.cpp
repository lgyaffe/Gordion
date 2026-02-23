#include "Print.h"
#include "Rep.h"
#include "Build.h"
#include "Canon.h"
#include "Numerics.h"
#include "Save.h"
#include "Blab.h"
#include <locale>
#include <cstdio>
#include <regex>

void Print::print_obs (uint i, uint j)		// Print ObsList::obs range
    {
    if (i < ObsList::obs.nobs() && j < ObsList::obs.nobs())
	{
	for (int k(i) ; k <= j ; ++k) ObsList::obs.print (cout, k) ;
	}
    else cout << "Invalid observable number " << i << " or " << j << "\n" ;
    }

void Print::print_obs (uint i)			// Print indexed Obs
    {
    if (i < ObsList::obs.nobs())
	{
	ObsList::obs.print (cout, i) ;
	}
    else gripe (format("Invalid observable number {}", i)) ;
    }

void Print::print_obs (string word)		// Print specified Obs
    {
    bool	found { false } ;
    SymbStr	s { to_symb(word) } ;
    uint	indx { ObsList::obs.find (s) } ;
    if (indx != UINT_MAX)
	{
	ObsList::obs.print (cout, indx) ;
	found = true ;
	}
    if (!found) gripe (format("Observable {} not known", word)) ;
    }

void Print::print_obs ()			// Print ObsList::obs
    {
    ObsList::obs.print (cout) ;
    }

void Print::print_obs_select (string word)	// Print order-selected Obs
    {
    std::regex	patt3 { "\\((\\d+),(\\d+),([A-Z]+[a-z]+)\\)" } ;
    std::regex	patt2 { "\\((\\d+),(\\d+)\\)" } ;
    std::regex	patt1 { "\\(([A-Z]+[a-z]+)\\)" } ;
    std::smatch match ;
    if (std::regex_match (word, match, patt3))
	{
	int	cord { stoi (match[1].str()) } ;
	int	xord { stoi (match[2].str()) } ;
	string	otyp { match[3].str() } ;
	auto ptr { Obs::obstypes.find (otyp) } ;
	if (ptr != Obs::obstypes.end())
	    {
	    ObsType type { ptr->second } ;
	    for (int k(0) ; k < ObsList::obs.nobs() ; ++k)
		{
		const Obs& obs { ObsList::obs(k) } ;
		if (obs.corder == cord && obs.xorder == xord && obs.type == type)
		    {
		    ObsList::obs.print (cout, k) ;
		    }
		}
	    }
	else gripe (format("Invalid obs type {}",otyp)) ;
	}
    else if (std::regex_match (word, match, patt2))
	{
	int	cord { stoi (match[1].str()) } ;
	int	xord { stoi (match[2].str()) } ;
	for (int k(0) ; k < ObsList::obs.nobs() ; ++k)
	    {
	    const Obs& obs { ObsList::obs(k) } ;
	    if (obs.corder == cord && obs.xorder == xord)
		{
		ObsList::obs.print (cout, k) ;
		}
	    }
	}
    else if (std::regex_match (word, match, patt1))
	{
	string	otyp { match[1].str() } ;
	auto ptr { Obs::obstypes.find (otyp) } ;
	if (ptr != Obs::obstypes.end())
	    {
	    ObsType type { ptr->second } ;
	    for (int k(0) ; k < ObsList::obs.nobs() ; ++k)
		{
		const Obs& obs { ObsList::obs(k) } ;
		if (obs.type == type)
		    {
		    ObsList::obs.print (cout, k) ;
		    }
		}
	    }
	else gripe (format("Invalid obs type {}",otyp)) ;
	}
    else gripe (format ("Don't understand {}",word)) ;
    }

void Print::print_base ()			// Print ObsList::base
    {
    ObsList::base.print (cout) ;
    }

void Print::print_op (uint i, uint j)		// Print Op::list range
    {
    if (i < Op::list.size() && j < Op::list.size())
	{
	for (int k(i) ; k <= j ; ++k) Op::print (cout, k) ;
	}
    else gripe (format("Invalid operator number {} or j", i, j)) ;
    }

void Print::print_op (uint i)			// Print indexed Op
    {
    if (i < Op::list.size())
	{
	Op::print (cout, i) ;
	}
    else gripe (format("Invalid operator number {}", i)) ;
    }

void Print::print_op ()				// Print Op::list
    {
    Op::print (cout) ;
    }

void Print::print_primary ()			// Print primary Op's
    {
    cout << "Primary Op's:\n" ;
    for (const auto& op : Op::list)
	{
	if (op.primary) cout << op << "\n" ;
	}
    }

void Print::print_fermiinit ()			// Print fermi init map
    {
    const auto& list { ObsList::obs } ;
    cout << "Fermion initializations:\n" ;
    for (const auto& [indx_f,indx_g] : ObsList::fermiinit)
	{
	cout << "<" << list(indx_f) << "> = -<"
		    << list(indx_g) << ">/2\n" ;
	}
    }

void Print::print_gen (uint i, bool full)	// Print specified Gen
    {
    const auto&	gens { global.info().gens[global.repnum] } ;
    if (i < gens.size())
	{
	const Gen& gen { gens[i] } ;

	cout << Rep::list[global.repnum].name << " generator "
	     << global.fg() << i << (gen.T_odd ? "*" : "")
	     << " (" << gen.order << ") = " ;
	if (full) cout << gen ;
	else	  gen.print (cout) ;
	cout << "\n" ;
	}
    else gripe (format("Invalid generator number {}", i)) ;
    }

void Print::print_gen (bool full)		// Print all Gens
    {
    const auto&	gens { global.info().gens[global.repnum] } ;

    for (int i(0) ; i < gens.size() ; ++i) print_gen (i,full) ;
    }

void Print::print_rep (const string& word)	// Print Rep(s)
    {
    if (word == "*")
	{
	cout << "Defined Symmetry Projectors:\n" ;
	for (const auto& rep : Rep::list) cout << rep ;
	}
    else
	{
	try {
	    uint i { Rep::known (word) } ;
	    cout << Rep::list[i] ;
	    }
	catch (const exception& e)
	    {
	    gripe ("Unknown representation: " + word) ;
	    }
	}
    }

void Print::print_rep ()			// Report active Rep
    {
    cout << "Defined Irreps:" ;
    for (const auto& rep : Rep::list) cout << " " << rep.name ;
    cout << "\nActive rep: " << Rep::list[global.repnum].name << "\n" ;
    }

void Print::print_symm ()			// Print Symm::list
    {
    cout << "Symmetries:" ;
    for (const auto& symm : Symm::list) { cout << "\n  " << symm ; }
    cout << "\n" ;
    }

void Print::print_symm (const string& word)	// Print specified Symm
    {
    try {
	cout << Symm::list[Symm::known(std::move(word)).item] << "\n" ;
	}
    catch (const exception& e)
	{
	gripe ("Unknown symetry: " + word) ;
	}
    }

void Print::print_symmsets ()			// Print Canon::symmset
    {
    for (int i(0) ; i < Canon::symmset.size() ; ++i)
	{
	cout << "symmset[" << i << "] =" ;
	for (auto k : Canon::symmset[i])
	    {
	    cout << " " << Symm::list[k].name ;
	    }
	cout << "\n" ;
	}
    }

void Print::print_couplings ()			// Print Coupling::couplings
    {
    cout << "Couplings:\n" ;
    for (const auto& c : Coupling::list)
	{
	cout << "  " << c << " = " << c.value << "\n" ;
	}
    }

void Print::print_hamiltonian ()		// Print Hamiltonian
    {
    cout << "  H_" << global.fg() << " = "
	 << global.info().Hterms << "\n" ;
    }

void Print::print_freeenergy ()			// Print free energy
    {
    cout << "  F_" << global.fg() << " = "
	 << global.info().Hterms << "\n" ;
    }

void Print::print_grad (uint i, uint j)		// Print gradient element
    {
    const auto& grad { global.data().grad } ;

    if (i < grad.entry().ncol)
	{
	const auto& term { global.info().Hterms[i].cpoly } ;

	if (j < grad.entry().nrow)
	    {
	    if (grad(j,i).len)
		{
		cout << "grad (" << i << "," << j << ") ="
		     << " d(" << term << ")"
		     << "/d(" << global.fg() << j << ") = "
		     << grad(j,i) << "\n" ;
		}
	    }
	else gripe (format("Invalid generator number {}", j)) ;
	}
    else gripe (format("Invalid term number {}", i)) ;
    }

void Print::print_grad (uint i)			// Print gradient row
    {
    const auto& grad { global.data().grad } ;

    if (i < grad.entry().ncol)
	{
	for (int j(0) ; j < grad.entry().nrow ; ++j) print_grad (i,j) ;
	}
    else gripe (format("Invalid term number {}", i)) ;
    }

void Print::print_grad ()			// Print gradient
    {
    for (const auto& poly : global.data().grad)
	{
	const GradHdr&	info ( poly ) ;
	auto		i    { info.Hterm } ;
	auto		j    { info.gen   } ;
	const auto&	term { global.info().Hterms[i].cpoly } ;

	if (info.len)
	    {
	    cout << "grad (" << i << "," << j << ") ="
		 << " d(" << term << ")/d(" << global.fg() << j << ") = "
		 << poly << "\n" ;
	    }
	}
    }

void Print::print_curv (uint i, uint j, uint k)	// Print curvature element
    {
    const auto& curv { global.data().curv[global.repnum] } ;

    if (i < curv.entry().nslice)
	{
	const auto& term { global.info().Hterms[i].cpoly } ;

	if (j < curv.entry().ncol && k < curv.entry().nrow)
	    {
	    if (curv(k,j,i).len)
		{
		cout << Rep::list[global.repnum].name
		     << " curv (" << i << "," << j << "," << k << ") ="
		     << " d^2(" << term << ")"
		     << "/d(" << global.fg() << j << ")"
		     << "/d(" << global.fg() << k << ") = "
		     << curv(k,j,i) << "\n" ;
		}
	    }
	else gripe (format("Invalid generator number {} or {}", j, k)) ;
	}
    else gripe (format("Invalid term number {}", i)) ;
    }

void Print::print_curv (uint i, uint j)		// Print curvature matrix row
    {
    const auto& curv { global.data().curv[global.repnum] } ;

    if (i < curv.entry().nslice)
	{
	if (j < curv.entry().ncol)
	    {
	    for (int k(0) ; k < curv.entry().nrow ; ++k) print_curv (i,j,k) ;
	    }
	else gripe (format("Invalid generator number {}", j)) ;
	}
    else gripe (format("Invalid term number {}", i)) ;
    }

void Print::print_curv (uint i)			// Print curvature slice
    {
    const auto& curv { global.data().curv[global.repnum] } ;

    if (i < curv.entry().nslice)
	{
	for (int j(0) ; j < curv.entry().ncol ; ++j) print_curv (i,j) ;
	}
    else gripe (format("Invalid term number {}", i)) ;
    }

void Print::print_curv ()			// Print curvature
    {
    for (const auto& poly : global.data().curv[global.repnum])
	{
	const CurvHdr&	info ( poly ) ;
	auto		i    { info.Hterm } ;
	auto		j    { info.gen1  } ;
	auto		k    { info.gen2  } ;
	const auto&	term { global.info().Hterms[i].cpoly } ;

	if (info.len)
	    {
	    cout << Rep::list[global.repnum].name
		 << " curv (" << i << "," << j << "," << k << ") ="
		 << " d^2(" << term << ")/d(" << global.fg() << j 
		 << ")/d(" << global.fg() << k << ") = " << poly << "\n" ;
	    }
	}
    }

void Print::print_lagrange (uint i, uint j)	// Print Lagrange mtx element
    {
    const auto&	lagr  { global.data().lagr[global.repnum] } ;

    if (i < lagr.entry().ncol && j < lagr.entry().nrow)
	{
	if (lagr(j,j).len)
	    {
	    cout << Rep::list[global.repnum].name << " lagrange (" 
		 << global.fg() << i << ","
		 << global.fg() << j << ") = "
		 << lagr(j,i) << "\n" ;
	    }
	}
    else gripe (format("Invalid generator number {} or {}", i,j)) ;
    }

void Print::print_lagrange (uint i)		// Print Lagrange matrix row
    {
    const auto&	lagr  { global.data().lagr[global.repnum] } ;

    if (i < lagr.entry().ncol)
	{
	for (int j(0) ; j < lagr.entry().nrow ; ++j) print_lagrange (i,j) ;
	}
    else gripe (format("Invalid generator number {}", i)) ;
    }

void Print::print_lagrange ()			// Print Lagrange bracket
    {
    for (const auto& poly : global.data().lagr[global.repnum])
	{
	const LagrHdr&	info ( poly ) ;
	auto		i    { info.gen1 } ;
	auto		j    { info.gen2 } ;

	if (info.len)
	    {
	    cout << Rep::list[global.repnum].name << " lagrange ("
		 << global.fg() << i << ","
		 << global.fg() << j << ") = " << poly << "\n" ;
	    }
	}
    }

void Print::print_mode (uint i)
    {
    const auto& repnam { Rep::list[numerics.lastrep].name } ;
    const auto& modes  { numerics.modes } ;

    if (i < numerics.modes.n_cols)
	{
	cout << repnam << " mode " << i << " =\n" << modes.col(i) ;
	}
    }

void Print::print_mode ()
    {
    for (int i(0) ; i < numerics.modes.n_cols ; ++i) print_mode (i) ;
    }

void Print::print_spectrum ()
    {
    const auto& repnam   { Rep::list[numerics.lastrep].name } ;
    auto	prevprec { cout.precision(12) } ;
    cout << repnam << " spectrum =\n" ;
    numerics.spectrum.raw_print (cout) ;
    cout << std::setprecision (prevprec) ;
    }

void Print::print_geodesic (uint i, uint j)	// Print specified geodesic equation
    {
    const auto&	list  { ObsList::obs } ;

    if (i < global.nobs())
	{
	auto [stage,bckt,pos]	{ global.bckt_pos (i) } ;
	const auto& geos	{ global.data(stage).geos[bckt] } ;

	if (j < geos.entry().nrow)
	    {
	    if (global.geoswap) Save::read_geo_bckt (bckt) ;
	    if (geos(j,pos).len)
		{
		cout << " geo (" << i << "," << j << ") = d<"
		     << *list[i] << ">/d("
		     << global.fg() << j << ") = "
		     << geos(j,pos) << "\n" ;
		}
	    if (global.geoswap) Save::read_geo_bckt (-bckt-1) ;
	    }
	else gripe (format("Invalid generator number {}", j)) ;
	}
    else gripe (format("Invalid observable number {}", i)) ;
    }

void Print::print_geodesic (uint i)		// Print geo eqns for specified Obs
    {
    if (i < global.nobs())
	{
	auto [stage,bckt,pos]	{ global.bckt_pos (i) } ;
	const auto& geos	{ global.data(stage).geos[bckt] } ;

	for (int j(0) ; j < geos.entry().nrow ; ++j) print_geodesic (i,j) ;
	}
    else gripe (format("Invalid observable number {}", i)) ;
    }

void Print::print_geodesic ()			// Print all geodesic equations
    {
    auto& bckt { global.info().bckt } ;

    for (int bcktnum(0) ; bcktnum < bckt.size() ; ++bcktnum)
	{
	auto& geos { global.data().geos[bcktnum] } ;
	if (global.geoswap) Save::read_geo_bckt (bcktnum) ;

	for (const auto& poly : geos)
	    {
	    const GeoHdr&	info ( poly ) ;
	    auto		i    { info.indx } ;
	    auto		j    { info.gen } ;

	    if (info.len)
		{
		cout << " geo (" << i << "," << j << ") = "
		     << "d<" << ObsList::obs(i) << ">/d("
		     << global.fg() << j << ") = " << poly << "\n" ;
		}
	    }
	if (global.geoswap) Save::read_geo_bckt (-bcktnum-1) ;
	}
    }

void Print::print_cache ()			// Print canonicalization cache
    {
    cout << Canon::cache ;
    }

void Print::print_sysindex ()			// Print sys data index
    {
    if (global.sysfile.size())
	{
	cout << "RecID  nslice  ncol      nrow     reclen      filebeg      fileend\n" ;
	for (auto& entry : global.sysindex)
	    {
	    if (!entry.reclen) continue ;
	    ulong fileend { entry.filepos + entry.reclen * sizeof (RecHdr) } ;
	    cout << RecIndx::idname[(int) entry.id] << "\t"
		 << std::setw(3)  << (int) entry.nslice << " "
		 << std::setw(7)  << entry.ncol << " "
		 << std::setw(9)  << entry.nrow << " "
		 << std::setw(10)  << entry.reclen << " "
		 << std::setw(12) << entry.filepos << " "
		 << std::setw(12) << fileend << "\n" ;
	    }
	}
    else gripe ("No opened sys info file") ;
    }

void Print::print_vevindex ()			// Print vev data index
    {

    if (global.vevstream.is_open())
	{
	Couplings list { Coupling::list.size() } ;
	for (int i(0) ; Save::read_coup (i,&list) ; ++i)
	    {
	    if (!i) cout << "Vev data sets in " << global.vevfile << ":\n" ;
	    cout << "  #" << i ;
	    auto sep { ": " } ;
	    for (auto& c : list)
		{
		cout << sep << c << " = " << c.value ;
		sep = ", " ;
		}
	    cout << "\n" ;
	    }
	}
    else gripe ("No open vev data file") ;
    }

void Print::print_theory()			// Print theory
    {
    cout << "Theory is " << theory.name.data() << std::endl ;
    if (theory.box.point)
	{
	cout << "Lattice size: " << theory.box_size() << "\n" ;
	}
    }

void Print::print_version()			// Print program version
    {
    cout << "Program version: " << global.version << "\n" ;
    }

void Print::print_rkmethods()			// Print Runge-Kutta method name
    {
    cout << "Avaliable RK integration methods:\n" ;
    for (const auto& rk : RKdef::list)
	{
	cout << "  " << rk.name << "\n" ;
	}
    cout << "Curren RK method: " << numerics.rk.name  << "\n" ;
    }

void Print::print_stats ()			// Print global statistics
    {
    cout.imbue (std::locale("en_US.UTF-8")) ;
    cout << "  Canonical obs:        " << ObsList::obs.nobs()       << "\n" ;
    cout << "  Basic obs:            " << ObsList::base.nobs()      << "\n" ;
    cout << "  Operator terms:       " << Op::list.size()           << "\n" ;
    cout << "  Ode integrations:     " << Ode::stats.integs         << "\n" ;
    cout << "      steps:            " << Ode::stats.steps          << "\n" ;
    cout << "      rejects:          " << Ode::stats.rejects        << "\n" ;
    cout << "      derivs:           " << Ode::stats.derivs         << "\n" ;
    cout << "      failures:         " << Ode::stats.fails          << "\n" ;
    cout << "  Minimizations:        " << numerics.stats.tries      << "\n" ;
    cout << "      failures:         " << numerics.stats.fails      << "\n" ;
    cout << "  Commutations:         " << global.count().commutes   << "\n" ;
    cout << "  Obs canonicalized:    " << global.count().canons     << "\n" ;
    cout << "      assessed:         " << global.count().assessed   << "\n" ;
    cout << "      found             " << global.count().found      << "\n" ;
    cout << "      bounded:          " << global.count().bounded    << "\n" ;
    cout << "      discarded:        " << global.count().discarded  << "\n" ;
    cout << "      classified:       " << global.count().classified << "\n" ;
    cout << "      factored:         " << global.count().factored   << "\n" ;
    cout << "      reduced:          " << global.count().reduced    << "\n" ;
    cout << "      stored:           " << global.count().stored     << "\n" ;
    cout << "      approximated:     " << global.count().approxed   << "\n" ;
    cout << "      cannot approx:    " << global.count().noapprox   << "\n" ;
    cout << "      no intersect:     " << global.count().nointersect<< "\n" ;
    cout << "  Symmetry transforms:  " << global.count().symmtrans  << "\n" ;
    cout << "           elements:    " << Symm::list.size()         << "\n" ;
    cout << "           subsets:     " << Canon::stats.symmsubsets  << "\n" ;
    cout << "           irreps:      " << Rep::list.size()          << "\n" ;
    cout << "  Canon cache hits:     " << global.count().cachehits  << "\n" ;
    cout << "              misses:   " << global.count().cachemiss  << "\n" ;
    cout << "              size:     " << Canon::cache.size()       << "\n" ;
    cout << "  Loop table size:      " << Canon::looptable.size()   << "\n" ;
    cout << "             nodes:     " << Canon::stats.looptblnodes << "\n" ;
    cout << "  Spec table size:      " << Canon::spectable.size()   << "\n" ;
    cout << "             nodes:     " << Canon::stats.spectblnodes << "\n" ;
    cout.imbue (std::locale::classic()) ;
    }

void Print::print_state ()			// Print global state variables
    {
    int	stage ( global.stage ) ;
    cout << " Theory:             " << theory.name.data() << "\n" ;
    cout << " Lattice size:       " << theory.box_size() << "\n" ;;
    cout << " Stage:              " << (stage ? "fermion" : "gauge") << "\n" ;
    cout << " Active rep:         " << Rep::list[global.repnum].name << "\n" ;
    cout << " Max obs order:      " << global.maxord()   << "\n" ;
    cout << " Max gen order:      " << global.maxgen()   << "\n" ;
    cout << " Obs approx:         " << global.approx     << "\n" ;
    cout << " Max threads:        " << global.maxthread  << "\n" ;
    cout << " Auto save:          " << global.autosave   << "\n" ;
    cout << " Neg curvature ok:   " << global.oknegeig   << "\n" ;
    cout << " Symmetrize curv:    " << global.symcurv    << "\n" ;
    cout << " Auto E-loop gens:   " << Gen::autoEgens    << "\n" ;
    cout << " Check obs:          " << Obs::check        << "\n" ;
    cout << " Dot obs:            " << SymbStr::dots     << "\n" ;
    cout << " Geo swap:           " << global.geoswap    << "\n" ;
    cout << " Normalize gens:     " << Gen::gennorm      << "\n" ;
    cout << " Max Newton iters:   " << numerics.minmax   << "\n" ;
    cout << " Max ODE steps:      " << numerics.odemax   << "\n" ;
    cout << " Minimize gen limit: " << numerics.minlim   << "\n" ;
    cout << " Minimize tolerance: " << numerics.mintol   << "\n" ;
    cout << " Ode tolerance:      " << numerics.odetol   << "\n" ;
    cout << " Ode RK method:      " << numerics.rk.name  << "\n" ;
    cout << " Svd truncation:     " << numerics.svdlim   << "\n" ;
    cout << " Spectrum gen limit: " << numerics.speclim  << "\n" ;
    cout << " Save directory:     " << global.savedir    << "\n" ;
    cout << " Sys info file:      " << global.sysfilename() << "\n" ;
    cout << " Vev data file:      " << global.vevfilename() << "\n" ;
    cout << " Vev file append:    " << global.vevappend     << "\n" ;
    cout << " MMA directory:      " << global.MMAdir        << "\n" ;
    cout << " MMA result file:    " << global.MMAfilename() << "\n" ;
    cout << " MMA file append:    " << global.MMAappend     << "\n" ;
    cout << " MMA save limit:     " << global.info().MMAlimit << "\n" ;
    cout << " MMA addl save list: " << global.info().MMAlist  << "\n" ;
    }

void Print::print_blab ()			// Print blab levels
    {
    for (const auto& [name,indx] : Blab::blabmap)
	{
	cout << " Blab level [" << name << "]\t= " << Blab::blablevel[indx] << "\n" ;
	}
    }

void Print::print_bcktlist ()			// Print bucket list
    {
    auto& bcktlist { global.info().bckt } ;

    cout << "    Bucket      Beg         End\n" ;
    for (auto& bckt : bcktlist)
	{
	cout << format(" {:6d} {:11d} {:11d}\n", bckt[0], bckt[1], bckt[2]) ;
	}
    }

void Print::print_obsstats ()			// Print observable statistics
    {
    ObsStats		obsstats { ObsList::obs } ;
    auto		maxc { obsstats.maxc } ;
    auto		maxx { obsstats.maxx } ;
    auto		maxl { obsstats.maxloop } ;
    ulong		total[]  { 0, 0, 0, 0, 0 } ;
    static ObsType	types[] { ObsType::Loop, ObsType::Eloop, ObsType::EEloop,
				  ObsType::Fermion, ObsType::Efermion } ;

    cout << " scorder corder xorder     Loop     E-loop    EE-loop    Fermion  E-fermion\n" ;

    for (int sc(0) ; sc <= maxc + maxx ; ++sc)
	{
	for (int c(0) ; c <= sc ; ++c)
	    {
	    int x { sc - c } ;
	    if (!obsstats.get(c,x)) continue ;

	    cout << format("{:6d} {:6d} {:6d}", sc, c, x) ;
	    for (auto t : types)
		{
		ulong count { obsstats.get((int) t,c,x) } ;
		if (count)	cout << format(" {:10d}", count) ;
		else		cout << format(" {:>10s}", "-") ;
		total[(int) t] += count ;
		}
	    cout << "\n" ;
	    }
	}
    cout << std::setw(20) << "" ;
    for (auto t : types)
	{
	ulong count { total[(int) t] } ;
	if (count) cout << format(" {: >10s}", "---") ;
	else	   cout << format(" {:10s}", "") ;
	}
    cout << format("\n {:19s}", "Totals by type:") ;
    for (auto t : types)
	{
	ulong count { total[(int) t] } ;
	if (count) cout << format(" {:10d}", count) ;
	else	   cout << format(" {:10s}", "") ;
	}
    cout << "\n" ;
    cout << format(" Avg Obs Length: {:<4.1f}\n", obsstats.avglen) ;
    if (maxl)
	{
	auto maxvev { numerics.vev[maxl] } ;
	cout << format(" Max Loop Vev: #{} = {:4.2f}\n", maxl, maxvev) ;
	}
    return ;
    }

void Print::print_geostats ()			// Print geodesic statistics
    {
    cout.imbue (std::locale("en_US.UTF-8")) ;
    cout << "Geodesic equations:\t" << global.count().ngeos << "\n" ;
    cout << "\t terms:   \t" << global.count().geoterms << "\n" ;
    for (int ord (0) ; ord <= PSIZ ; ++ord)
	{
	if (ulong k { global.count().geotermord [ord] })
	    cout << "\t   order " << ord << ":\t" << k << "\n" ;
	}
    cout.imbue (std::locale::classic()) ;
    }
