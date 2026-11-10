#include "Print.h"
#include "Canon.h"
#include "Gen.h"
#include "Numerics.h"
#include "Rep.h"
#include "Save.h"
#include "Blab.h"
#include <locale>
#include <cstdio>
#include <regex>

ostream& Print::coeffprt (ostream& stream, doub c)	// Pretty print coefficient
    {
    if (c)
	{
	const string sgn { c > 0 ? "+" : "-" } ;
	c = std::abs(c) ;
	doub x { c * sqrt(2.0) } ;
	if      (is_one (c))	stream << sgn ;
	else if (is_int (c))	stream << sgn <<    c << " " ;
	else if (is_int (1/c))	stream << sgn << "1/" << (int)(1/c) << " " ;
	else if (is_int (2*c))	stream << sgn << 2*c  << "/2" << " " ;
	else if (is_int (3*c))	stream << sgn << 3*c  << "/3" << " " ;
	else if (is_int (4*c))	stream << sgn << 4*c  << "/4" << " " ;
	else if (is_int (6*c))	stream << sgn << 6*c  << "/6" << " " ;
	else if (is_int (8*c))	stream << sgn << 8*c  << "/8" << " " ;
	else if (is_int (9*c))	stream << sgn << 9*c  << "/9" << " " ;
	else if (is_int (12*c))	stream << sgn << 12*c << "/12" << " " ;
	else if (is_int (16*c))	stream << sgn << 16*c << "/16" << " " ;
	else if (is_int (x))	stream << sgn <<   x  << "/\u221A2" << " "  ;
	else if (is_int (2*x))	stream << sgn << 2*x  << "/(2\u221A2)" << " " ;
	else if (is_int (3*x))	stream << sgn << 3*x  << "/(3\u221A2)" << " " ;
	else if (is_int (4*x))	stream << sgn << 4*x  << "/(4\u221A2)" << " " ;
	else if (is_int (8*x))	stream << sgn << 8*x  << "/(8\u221A2)" << " " ;
	else if (is_int (16*x))	stream << sgn << 16*x << "/(16\u221A2)" << " " ;
	else			stream << sgn << format("{:.2e} ", c) ;
	}
    else stream << "+0 " ;
    return stream ;
    } ;

void Print::print_obs (numb i, numb j)		// Print global.obs range
    {
    if (global.obs.swapped) Save::reload_obs() ;

    if (i < global.obs.size() && j < global.obs.size())
	{
	for (int k(i) ; k <= j ; ++k) global.obs.print (cout, k) ;
	}
    else cout << "Invalid observable number " << i << " or " << j << "\n" ;
    }

void Print::print_obs (numb i)			// Print indexed Obs
    {
    if (global.obs.swapped) Save::reload_obs() ;

    if (i < global.obs.size())
	{
	global.obs.print (cout, i) ;
	}
    else gripe (format("Invalid observable number {}", i)) ;
    }

void Print::print_obs (const string& word)	// Print specified Obs
    {
    if (global.obs.swapped) Save::reload_obs() ;
    numb indx { global.obs.find (Str (word)) } ;
    if (indx < MAXNUM) global.obs.print (cout, indx) ;
    else gripe (format("Observable {} not known", word)) ;
    }

void Print::print_obs ()			// Print global.obs
    {
    if (global.obs.swapped) Save::reload_obs() ;
    global.obs.print (cout) ;
    }

void Print::print_obs_select (const string& word)	// Print type-selected Obs
    {
    if (global.obs.swapped) Save::reload_obs() ;

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
	    for (int k(0) ; k < global.obs.size() ; ++k)
		{
		const Obs& obs { global.obs(k) } ;
		if (obs.corder == cord && obs.xorder == xord && obs.type == type)
		    {
		    global.obs.print (cout, k) ;
		    }
		}
	    }
	else gripe ("Invalid obs type " + otyp) ;
	}
    else if (std::regex_match (word, match, patt2))
	{
	int	cord { stoi (match[1].str()) } ;
	int	xord { stoi (match[2].str()) } ;
	for (int k(0) ; k < global.obs.size() ; ++k)
	    {
	    const Obs& obs { global.obs(k) } ;
	    if (obs.corder == cord && obs.xorder == xord)
		{
		global.obs.print (cout, k) ;
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
	    for (int k(0) ; k < global.obs.size() ; ++k)
		{
		const Obs& obs { global.obs(k) } ;
		if (obs.type == type)
		    {
		    global.obs.print (cout, k) ;
		    }
		}
	    }
	else gripe (format("Invalid obs type {}", otyp)) ;
	}
    else gripe ("Don't understand " + word) ;
    }

void Print::print_op (numb i, numb j)		// Print OpList range
    {
    const auto& oplist { global.info().ops } ;
    if (i < oplist.size() && j < oplist.size())
	{
	for (int k(i) ; k <= j ; ++k) oplist.print (cout, k) ;
	}
    else gripe (format("Invalid operator number {} or j", i, j)) ;
    }

void Print::print_op (numb i)			// Print indexed Op
    {
    const auto& oplist { global.info().ops } ;
    if (i < oplist.size()) oplist.print (cout, i) ;
    else gripe (format("Invalid operator number {}", i)) ;
    }

void Print::print_op ()				// Print OpList
    {
    const auto& oplist { global.info().ops } ;
    oplist.print (cout) ;
    }

void Print::print_gen (uint i, bool full)	// Print specified Gen
    {
    const auto&	gens { global.info().gens[global.repnum] } ;
    if (i < gens.size())
	{
	const Gen& gen { gens[i] } ;

	cout << "  " << Rep::list[global.repnum].name
	     << " generator " << global.fg() << i
	     << (gen.T_odd ? "*" : "") << " (" << gen.order << ") = " ;
	if (full) cout << gen ;
	else	  gen.print (cout) ;
	cout << "\n" ;
	}
    else gripe (format("Invalid generator number {}", i)) ;
    }

void Print::print_gen (bool full)		// Print all Gens
    {
    const auto&	gens { global.info().gens[global.repnum] } ;
    int		live (0) ;

    for (int i(0) ; i < gens.size() ; ++i)
	{
	if (!gens[i].active) continue ;
	if (!live++) cout << "Active generators:\n" ;
	print_gen (i,full) ;
	}
    if (live == gens.size()) return ;
    cout << "Inactive generators:\n" ;
    for (int i(0) ; i < gens.size() ; ++i)
	{
	if (!gens[i].active) print_gen (i,full) ;
	}
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

void Print::print_symm (bool detail)		// Print Symm::list
    {
    cout << "Symmetries:" ;
    for (const auto& symm : Symm::list)
	{
	if (detail)	cout << "\n  " << symm ;
	else		cout << " " << symm.name ;
	}
    cout << "\n" ;
    }

void Print::print_symm (const string& word)	// Print specified Symm
    {
    try {
	cout << Symm::list[Symm::knownindx(std::move(word))] << "\n" ;
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

void Print::print_theory()			// Print theory
    {
    cout << "Theory is " << theory.name.data() << std::endl ;
    if (theory.box.point)
	{
	cout << "Lattice size: " << theory.box_size() << "\n" ;
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

void Print::print_ham_or_free (int i)		// Print H/F element
    {
    const auto&	ham	{ global.data().ham } ;
    const auto& Hterms	{ global.info().Hterms } ;
    auto	nterms	{ ham.entry().nelem } ;
    string	name	{ theory.euclid ? "  F_" : "  H_" } ;
    int		end	( i < 0 ? nterms : i + 1 ) ;
    int		indx	{ i < 0 ? 0 : i } ;

    if (end > nterms) gripe (format("Invalid term number {}", i)) ;
    for (; indx < end ; ++indx)
	{
	if (ham(indx).len())
	    {
	    cout << name << global.fg() << " (" << indx << ") ="
		 << ham(indx) << "\n" ;
	    }
	}
    }

void Print::print_ham_or_free ()		// Print H or F
    {
    const auto&	ham	{ global.data().ham } ;
    const auto&	Hterms	{ global.info().Hterms } ;
    auto	nterms	{ ham.entry().nelem } ;
    auto	recptr	{ ham.begin() } ;
    string	name	{ theory.euclid ? "  F_" : "  H_" } ;
    string	sep ;

    cout << name << global.fg() << " = " << Hterms ;
    if (nterms)
	{
	cout << "\n      = " ;
	for (int i(0) ; i < nterms ; ++i, ++recptr)
	    {
	    const PolyElem* ptr { &(*recptr) } ;
	    cout << sep << Hterms[i].coeff << "(" << *ptr << " ) " ;
	    sep = "+ " ;
	    }
	}
    cout << "\n" ;
    }

void Print::print_grad (uint i, uint j)		// Print gradient element
    {
    const auto&	ham	{ global.data().ham  } ;
    const auto& grad	{ global.data().grad } ;
    auto	nterms	{ grad.entry().ncol } ;
    auto	neven	{ grad.entry().nrow } ;

    if (i < nterms)
	{
	if (j < neven)
	    {
	    if (grad(j,i).len())
		{
		cout << "grad (" << i << "," << j << ") ="
		     << " d(" << ham(i) << ")"
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
    const auto& grad	{ global.data().grad } ;
    auto	nterms	{ grad.entry().ncol } ;
    auto	neven	{ grad.entry().nrow } ;

    if (i < nterms)
	{
	for (int j(0) ; j < neven ; ++j) print_grad (i,j) ;
	}
    else gripe (format("Invalid term number {}", i)) ;
    }

void Print::print_grad ()			// Print gradient
    {
    const auto&	ham	{ global.data().ham  } ;
    const auto&	grad	{ global.data().grad } ;
    auto	nterms	{ grad.entry().ncol } ;
    auto	neven	{ grad.entry().nrow } ;
    auto	recptr	{ grad.begin() } ;

    for (int i(0) ; i < nterms ; ++i)
	{
	for (int j(0) ; j < neven ; ++j, ++recptr)
	    {
	    const PolyElem* ptr { &(*recptr) } ;
	    if (ptr->len())
		{
		cout << "grad (" << i << "," << j << ") ="
		     << " d(" << ham(i) << ")/d(" << global.fg() << j << ") = "
		     << *ptr << "\n" ;
		}
	    }
	}
    }

void Print::print_curv (uint i, uint j, uint k)	// Print curvature element
    {
    const auto&	ham	{ global.data().ham } ;
    const auto& curv	{ global.data().curv[global.repnum] } ;
    auto	nterms	{ curv.entry().nslice } ;
    auto	neven	{ curv.entry().ncol } ;

    if (curv.entry().nrow != neven) fatal ("Corrupted curv record!") ;

    if (i < nterms)
	{
	if (j < neven && k < neven)
	    {
	    if (curv(k,j,i).len())
		{
		cout << Rep::list[global.repnum].name
		     << " curv (" << i << "," << j << "," << k << ") ="
		     << " d^2(" << ham(i) << ")"
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
    const auto& curv	{ global.data().curv[global.repnum] } ;
    auto	nterms	{ curv.entry().nslice } ;
    auto	neven	{ curv.entry().ncol } ;

    if (i < nterms)
	{
	if (j < neven)
	    {
	    for (int k(0) ; k < neven ; ++k) print_curv (i,j,k) ;
	    }
	else gripe (format("Invalid generator number {}", j)) ;
	}
    else gripe (format("Invalid term number {}", i)) ;
    }

void Print::print_curv (uint i)			// Print curvature slice
    {
    const auto& curv	{ global.data().curv[global.repnum] } ;
    auto	nterms	{ curv.entry().nslice } ;
    auto	neven	{ curv.entry().ncol } ;

    if (i < nterms)
	{
	for (int j(0) ; j < neven ; ++j) print_curv (i,j) ;
	}
    else gripe (format("Invalid term number {}", i)) ;
    }

void Print::print_curv ()			// Print curvature
    {
    const auto&	ham	{ global.data().ham } ;
    const auto&	curv	{ global.data().curv[global.repnum] } ;
    auto	nterms	{ curv.entry().nslice } ;
    auto	neven	{ curv.entry().ncol } ;
    auto	recptr	{ curv.begin() } ;

    if (curv.entry().nrow != neven) fatal ("Corrupted curv record!") ;

    for (int i(0) ; i < nterms ; ++i)
	{
	for (int j(0) ; j < neven ; ++j)
	    {
	    for (int k(0) ; k < neven ; ++k, ++recptr)
		{
		const PolyElem* ptr { &(*recptr) } ;
		if (ptr->len())
		    {
		    cout << Rep::list[global.repnum].name
			 << " curv (" << i << "," << j << "," << k << ") ="
			 << " d^2(" << ham(i) << ")/d(" << global.fg() << j 
			 << ")/d(" << global.fg() << k << ") = " << *ptr << "\n" ;
		    }
		}
	    }
	}
    }

void Print::print_lagrange (uint j, uint k)	// Print Lagrange mtx element
    {
    const auto&	lagr	{ global.data().lagr[global.repnum] } ;
    auto	neven	{ lagr.entry().ncol } ;
    auto	nodd	{ lagr.entry().nrow } ;

    if (j < neven && k-neven >= 0 && k-neven < nodd)
	{
	if (lagr(k-neven,j).len())
	    {
	    cout << Rep::list[global.repnum].name << " lagrange (" 
		 << global.fg() << j << ","
		 << global.fg() << k << ") = "
		 << lagr(k-neven,j) << "\n" ;
	    }
	}
    else gripe (format("Invalid generator number {} or {}", j,k)) ;
    }

void Print::print_lagrange (uint j)		// Print Lagrange matrix row
    {
    const auto&	lagr	{ global.data().lagr[global.repnum] } ;
    auto	neven	{ lagr.entry().ncol } ;
    auto	nodd	{ lagr.entry().nrow } ;

    if (j < neven)
	{
	for (int k(0) ; k < nodd ; ++k) print_lagrange (j,k+neven) ;
	}
    else gripe (format("Invalid generator number {}", j)) ;
    }

void Print::print_lagrange ()			// Print Lagrange bracket
    {
    const auto&	lagr	{ global.data().lagr[global.repnum] } ;
    auto	neven	{ lagr.entry().ncol } ;
    auto	nodd	{ lagr.entry().nrow } ;
    auto	recptr	{ lagr.begin() } ;

    for (int j(0) ; j < neven ; ++j)
	{
	for (int k(0) ; k < nodd ; ++k, ++recptr)
	    {
	    const PolyElem* ptr { &(*recptr) } ;
	    if (ptr->len())
		{
		cout << Rep::list[global.repnum].name << " lagrange ("
		     << global.fg() << j << ","
		     << global.fg() << k+neven << ") = " << *ptr << "\n" ;
		}
	    }
	}
    }

void Print::print_spectrum ()			// Print spectrum eigenvalues
    {
    auto prevprec { cout.precision(12) } ;
    cout << Rep::list[numerics.lastrep].name ;
    numerics.spectrum.raw_print (" spectrum =") ;
    cout << setprecision (prevprec) ;
    }

void Print::print_mode (uint i)			// Print spectrum eigenvectors
    {
    const auto& repnam { Rep::list[numerics.lastrep].name } ;
    const auto& modes  { numerics.modes } ;

    if (i < numerics.modes.n_cols)
	{
	cout << repnam << " mode " << i << " =\n" << modes.col(i) ;
	}
    }

void Print::print_mode ()			// Print spectrum eigenvectors
    {
    for (int i(0) ; i < numerics.modes.n_cols ; ++i) print_mode (i) ;
    }

void Print::print_geodesic (numb i, uint j)	// Print specified geodesics
    {
    const auto&	list  { global.obs } ;

    if (list.swapped) Save::reload_obs() ;
    if (i < list.size())
	{
	auto [stage,bckt,pos]	{ global.bckt_pos (i) } ;
	const auto& geos	{ global.data(stage).geos[bckt] } ;

	if (j < geos.entry().nrow)
	    {
	    if (global.geoswap) Save::read_geo_bckt (stage, bckt) ;
	    if (geos(j,pos).len())
		{
		cout << " geo (" << i << "," << j << ") = d<"
		     << list(i) << ">/d("
		     << global.fg(stage) << j << ") = "
		     << geos(j,pos) << "\n" ;
		}
	    if (global.geoswap) Save::read_geo_bckt (stage, -bckt-1) ;
	    }
	else gripe (format("Invalid generator number {}", j)) ;
	}
    else gripe (format("Invalid observable number {}", i)) ;
    }

void Print::print_geodesic (numb i)		// Print geodesic for specified Obs
    {
    if (global.obs.swapped) Save::reload_obs() ;
    if (i < global.obs.size())
	{
	auto [stage,bckt,pos]	{ global.bckt_pos (i) } ;
	const auto& geos	{ global.data(stage).geos[bckt] } ;

	for (int j(0) ; j < geos.entry().nrow ; ++j) print_geodesic (i,j) ;
	}
    else gripe (format("Invalid observable number {}", i)) ;
    }

void Print::print_geodesic ()			// Print all geodesic equations
    {
    if (global.obs.swapped) Save::reload_obs() ;

    for (const auto& bckt : global.info().bckt)
	{
	const auto&	list	{ global.obs } ;
	const auto&	stage	{ global.stage } ;
	const auto&	bcktnum	{ bckt[0] } ;
	const auto&	first	{ bckt[1] } ;
	const auto&	last	{ bckt[2] } ;
	const auto&	geos	{ global.data().geos[bcktnum] } ;
	auto		ngeos	{ geos.entry().ncol } ;
	auto		neven	{ geos.entry().nrow } ;

	if (global.geoswap) Save::read_geo_bckt (stage, bcktnum) ;

	auto		recptr	{ geos.begin() } ;
	for (int i(first) ; i <= last ; ++i)
	    {
	    for (int j(0) ; j < neven ; ++j, ++recptr)
		{
		const PolyElem* ptr { &(*recptr) } ;
		if (ptr->len())
		    {
		    cout << " geo (" << i << "," << j << ") = " << "d<" << list(i)
			 << ">/d(" << global.fg() << j << ") = " << *ptr << "\n" ;
		    }
		}
	    }
	if (global.geoswap) Save::read_geo_bckt (stage, -bcktnum-1) ;
	}
    }

void Print::print_base ()			// Print basic Obs
    {
    global.baseobs.print (cout) ;
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

void Print::print_blab ()			// Print blab levels
    {
    for (const auto& [name,indx] : Blab::blabmap)
	{
	cout << " Blab level [" << name << "]\t= " << Blab::level(indx) << "\n" ;
	}
    }

void Print::print_canon ()			// Print canonicalization cache
    {
    using Strpair = pair<Str,int> ;
    const auto& cache { Canon::cache } ;
    auto mysort { [](const Strpair a, const Strpair b)
		{ return Str::compare (a.first, b.first) ; } } ;
    cout << "Short Obs cache:\n" ;
    vector<Strpair> items ;
    items.reserve (cache.size()) ;
    for (const auto& p : cache) items.push_back (p) ;
    std::sort (items.begin(), items.end(), mysort) ;
    for (const auto& [s,indx] : items)
	{
	cout << " " << std::left << std::setw(8) << Str(s).print()
	     << " -> " << std::right << (indx < 0 ? '-' : ' ')
	     << global.obs (std::abs(indx)) << "\n" ;
	}
    }

void Print::print_fermiinit ()			// Print fermi init map
    {
    if (global.obs.swapped) Save::reload_obs() ;

    const auto& list { global.obs } ;
    cout << "Fermion initializations:\n" ;
    for (const auto& [indx_f,indx_g] : global.obs.fermiinit)
	{
	cout << "<" << list(indx_f) << "> = -<"
		    << list(indx_g) << ">/2\n" ;
	}
    }

void Print::print_obsstats ()			// Print observable statistics
    {
    if (global.obs.swapped) Save::reload_obs() ;

    const auto&		types	 { Obs::type_name } ;
    ObsStats		obsstats { global.obs } ;
    auto&		maxl	 { obsstats.maxloop } ;
    vector<ulong>	total ;

    cout << " scorder corder xorder" ;
    for (int t(0) ; t < obsstats.size() ; ++t)
	{
	total.push_back (obsstats.get(t)) ;
	if (total.back()) cout << std::setw(11) << types[t] ;
	}
    cout << "\n" ;

    for (int sc(0) ; sc <= obsstats.maxc + obsstats.maxx ; ++sc)
	{
	for (int c(0) ; c <= sc ; ++c)
	    {
	    int x { sc - c } ;
	    if (!obsstats.get(c,x)) continue ;
	    cout << format("{:6d} {:6d} {:6d} ", sc, c, x) ;
	    for (int t(0) ; t < obsstats.size() ; ++t)
		{
		if (!total[t]) continue ;
		auto count { obsstats.get(t,c,x) } ;
		if (count)	cout << std::setw(11) << count ;
		else		cout << std::setw(11) << "-" ;
		}
	    cout << "\n" ;
	    }
	}
    cout << std::setw(21) << "" ;
    for (auto& count : total) if (count) cout << std::setw(11) << "---" ;
    cout << "\n" ;
    cout << std::left << std::setw(21) << " Totals by type:" << std::right ;
    for (auto& count : total) if (count) cout << std::setw(11) << count ;
    cout << format("\n Avg Obs Length: {:<4.1f}\n", obsstats.avglen) ;
    if (maxl)
	{
	auto maxvev { numerics.vev[maxl] } ;
	cout << format(" Max Loop Vev: #{} = {:4.2f}\n", maxl, maxvev) ;
	}
    return ;
    }

void Print::print_polystats ()			// Print polynomial statistics
    {
    try { cout.imbue (std::locale("en_US.UTF-8")) ; }
    catch (const std::exception&) {}
    const auto& gradstats { global.count().gradstats } ;
    const auto& curvstats { global.count().curvstats } ;
    const auto& lagrstats { global.count().lagrstats } ;
    const auto& geostats  { global.count().geostats  } ;

    cout << std::setw(16) << std::left << "Poly statistics:" << std::right
	 << std::setw(PW) << "eqns"
	 << std::setw(PW) << "terms"
	 << std::setw(PW-1) << "order 0" ;
    for (int i(1) ; i <= PSIZ ; ++i) cout << std::setw(PW-1) << "order " << i ;
    cout << "\n" ;
    if (gradstats.eqns()) cout << std::setw(16) << "Gradient:"  << gradstats ;
    if (curvstats.eqns()) cout << std::setw(16) << "Curvature:" << curvstats ;
    if (lagrstats.eqns()) cout << std::setw(16) << "Lagrange brkt:" << lagrstats ;
    if  (geostats.eqns()) cout << std::setw(16) << "Geodesic eqns:" << geostats  ;
    cout.imbue (std::locale::classic()) ;
    }

void Print::print_primary ()			// Print primary Op's
    {
    const auto& oplist { global.info().ops } ;
    cout << "Primary Op's:\n" ;
    for (const auto& op : oplist)
	{
	if (op.primary) cout << op << "\n" ;
	}
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

void Print::print_stage ()			// Print minimization stage
    {
    auto stage { global.stage ? "fermion" : "gauge" } ;
    cout << "Current minimization stage: " << stage << "\n" ;
    }

void Print::print_state ()			// Print global state variables
    {
    int	stage ( global.stage ) ;
    cout << " Theory:             " << theory.name.data() << "\n" ;
    cout << " Lattice size:       " << theory.box_size() << "\n" ;;
    cout << " Stage:              " << (stage ? "fermion" : "gauge") << "\n" ;
    cout << " Active rep:         " << Rep::list[global.repnum].name << "\n" ;
    cout << " Max obs order:      " << global.info().maxord << "\n" ;
    cout << " Max gen order:      " << global.info().maxgen << "\n" ;
    cout << " Max threads:        " << global.maxthread  << "\n" ;
    cout << " Obs approx:         " << global.approx     << "\n" ;
    cout << " Add curvature obs:  " << global.xtraobs    << "\n" ;
    cout << " Auto E-loop gens:   " << Gen::autoEgens    << "\n" ;
    cout << " Auto save:          " << Save::autosave    << "\n" ;
    cout << " Geo swap:           " << global.geoswap    << "\n" ;
    cout << " Obs swap:           " << global.obsswap    << "\n" ;
    cout << " Check obs:          " << Obs::check        << "\n" ;
    cout << " Dot obs:            " << Str::dots         << "\n" ;
    cout << " Normalize gens:     " << Gen::gennorm      << "\n" ;
    cout << " Max Newton iters:   " << numerics.maxnewt  << "\n" ;
    cout << " Max ODE steps:      " << numerics.maxode   << "\n" ;
    cout << " Ode tolerance:      " << numerics.odetol   << "\n" ;
    cout << " Minimize tolerance: " << numerics.mintol   << "\n" ;
    cout << " Neg curvature ok:   " << numerics.oknegeig << "\n" ;
    cout << " Symmetrize curv:    " << numerics.symcurv  << "\n" ;
    cout << " Ode RK method:      " << numerics.rk.name  << "\n" ;
    cout << " SVD cutoff:         " << numerics.svdcut   << "\n" ;
    cout << " Input directory:    " << global.inputdir   << "\n" ;
    cout << " Save directory:     " << global.savedir    << "\n" ;
    cout << " Sys info file:      " << global.info().sysfile.path   << "\n" ;
    cout << " Vev data file:      " << global.info().vevfile.path   << "\n" ;
    cout << " Vev file append:    " << global.info().vevfile.append << "\n" ;
    cout << " MMA directory:      " << global.MMAdir                << "\n" ;
    cout << " MMA result file:    " << global.info().MMAfile.path   << "\n" ;
    cout << " MMA file append:    " << global.info().MMAfile.append << "\n" ;
    cout << " MMA obs subset:     " << global.info().MMAfile.obs    << "\n" ;
    }

void Print::print_stats ()			// Print global statistics
    {
    try { cout.imbue (std::locale("en_US.UTF-8")) ; }
    catch (const std::exception&) {}
    cout << "  Basic obs:            " << global.baseobs.size()     << "\n" ;
    cout << "  Canonical obs:        " << global.obs.size()         << "\n" ;
    cout << "      # gauge obs:      " << global.obs.nobsG          << "\n" ;
    cout << "      # fermi obs:      " << global.obs.nobsF          << "\n" ;
    cout << "      gauge obs hash:   " << global.obs.hashG	    << "\n" ;
    cout << "      fermi obs hash:   " << global.obs.hashF	    << "\n" ;
    cout << "  Operators:            " << global.info(0).ops.size()
    					+ global.info(1).ops.size() << "\n" ;
    cout << "      gauge:            " << global.info(0).ops.size() << "\n" ;
    cout << "      fermion:          " << global.info(1).ops.size() << "\n" ;
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

void Print::print_sysindex ()			// Print sys data index
    {
    auto& info { global.info() } ;
    if (info.sysfile.path.size())
	{
	cout << "Save file: " << info.sysfile.path << "\n" ;
	cout << "RecID  nslice  ncol      nrow     reclen      filebeg      fileend\n" ;
	for (auto& entry : info.sysindex)
	    {
	    if (!entry.reclen) continue ;
	    auto fileend { entry.filepos + entry.reclen * sizeof (Element) } ;
	    cout << RecIndx::idname[(int) entry.id] << "\t"
		 << std::setw(3)  << (int) entry.nslice << " "
		 << std::setw(7)  << entry.ncol << " "
		 << std::setw(9)  << entry.nrow << " "
		 << std::setw(10) << entry.reclen << " "
		 << std::setw(12) << entry.filepos << " "
		 << std::setw(12) << fileend << "\n" ;
	    }
	}
    else gripe ("No opened sys info file") ;
    }

void Print::print_version()			// Print program version
    {
    cout << "Program version: " << global.version << "\n" ;
    }

void Print::print_vevindex ()			// Print vev data index
    {
    auto& info { global.info() } ;
    if (info.vevfile.stream.is_open())
	{
	int		i (0) ;
	Couplings*	listptr ;
	cout << "Vev data sets in " << info.vevfile.path << ":\n" ;
	while ((listptr = Save::read_coup (i,false)))
	    {
	    cout << "  #" << i++ ;
	    auto sep { ": " } ;
	    for (auto& c : *listptr)
		{
		cout << sep << c << " = " << c.value ;
		sep = ", " ;
		}
	    cout << "\n" ;
	    }
	}
    else gripe ("No open vev data file") ;
    }

