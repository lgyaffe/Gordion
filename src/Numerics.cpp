#include "Numerics.h"
#include "Data.h"
#include "Gen.h"
#include "Global.h"
#include "Print.h"
#include "Rep.h"
#include "Save.h"
#include "Blab.h"
#include "Gripe.h"
#include <filesystem>

void Numerics::do_flow (uint indx, doub v0, doub v1, doub inc)	// Flow coupling
    {
    auto&	coupling { Coupling::list[indx] } ;
    doub&	value    { coupling.value } ;
    int  	stage    { coupling.stage } ; 
    auto	prevprec { cout.precision(12) } ;

    if (!inc) gripe ("Must have non-zero coupling increment!") ;
    if ((v1 - v0)/inc < 0) inc *= -1 ;
    if (stage != global.stage) global.stageinit (stage) ;
    if (!stage && global.fermivev) initialize (1) ;

    cout << std::scientific ;
    for (value = v0 ;;)
	{
	cout << coupling << format(" = {:6.3f}:", value) << flush ;
	try {
	    do_minimize () ;
	    if (global.interrupt) return ;
	    write_data () ;
	    }
	catch (const Abort& e)
	    {
	    cout.precision (prevprec) ;
	    cout << std::defaultfloat << e.what() << "\n" ;
	    abort (format ("Aborting flow at {} = {:.6f}", coupling.data(), value)) ;
	    }
	if (value == v1)				 break ;
	else if (std::abs(value-v1) > std::abs(1.6*inc)) value += inc ;
	else if (std::abs(value-v1) > std::abs(inc))	 value  = (value + v1)/2 ;
	else						 value  = v1 ;
	}
    cout << std::setprecision (prevprec) << std::defaultfloat ;
    }

void Numerics::do_minimize ()			// Do minimization
    {
    int  steps ;
    doub tol   (0) ;
    uint iters (0) ;
    uint total (0) ;
    ++stats.tries ;
    status.reset() ;

    do  {
	if (global.interrupt || ++iters > maxnewt) break ;
	steps = do_step (tol) ;
	total += std::abs(steps) ;
	tol = mintol ;
	} while (steps > 0) ;

    status_rpt (iters,total) ;
    if (iters > maxnewt)
	{
	++stats.fails ;
	abort ("Minimization failed, max iterations exceeded") ;
	}
    else if (steps < 0)
	{
	++stats.fails ;
	abort ("ODE integration failed!") ;
	}
    if (global.autosave && !global.interrupt) Save::save_vev () ;
    }

int Numerics::do_step (doub tol)		// Do geodesic integration step
    {
    if (global.interrupt) return 0 ;

    const auto&	blab	{ Blab::level(Blab::NUMERICS) } ;
    Ode		ode	{ do_dvev, err_norm, odetol, rk, maxode } ;
    doub	dnorm	{ eval_delta() } ;
    long	n	{ global.info().nobs } ;
    doub	s (0) ;
    bool	ok ;

    if (blab > 1) cout << "  |delta| = " << dnorm << ", " ;
    if (blab > 2) delta.raw_print ("\n   delta =") ;

    if (tol && dnorm <= tol) return 0 ;			// converged?

    if (n != vev.size())				// integrate vev subvec
	{
	real*	vevptr { &vev [global.stage ? global.info(0).nobs : 0] } ; 
	Rvec	subvev ( vevptr, n, false, true ) ;
	ok = ode.integrate (s, 1.0, subvev) ;
	}
    else ok = ode.integrate (s, 1.0, vev) ;

    if (global.stage == Global::Fermi) global.fermivev = true ;
    else if (global.fermivev) initialize (1) ;

    if (blab > 1) cout << ode.steps << " step(s) "
		       << ode.rejects << " rejects "
		       << (theory.euclid ? "F = " : "H = ")
		       << eval_ham() << "\n" << flush ;

    return ok ? ode.steps : -ode.steps ;
    }

const Uvec& Numerics::eval_inuse (uint repnum, bool T_odd)
    {
    auto&	inuse	{ T_odd ? Todds : Tevens } ;
    const auto& gens	{ global.info().gens[repnum] } ;
    int		neven	{ global.info().neven[repnum] } ;
    int		ngen	( gens.size() ) ;
    int		beg	{ T_odd ? neven : 0 } ;
    int		end	{ T_odd ? ngen : neven } ;
    int		n(0) ;

    inuse.set_size(end-beg) ;
    for (int i(beg) ; i < end ; ++i)
	if (gens[i].active) inuse (n++) = i - beg ;
    inuse.resize(n) ;
    return inuse ;
    }

doub Numerics::eval_ham (bool print)		// Evaluate Hamiltonian/free energy
    {
    if (!global.maxord()) gripe ("Make observables first!") ;

    H = 0 ;
    for (const auto& Hterm : global.info().Hterms)
	{
	doub val (0.0) ;
	for (const auto& term : Hterm.cpoly) val += termvalue (term) ;
	H += Hterm.coeff() * val ;
	}
    if (print)
	{
	auto label { theory.euclid ? "Free energy" : "Hamiltonian" } ;
	auto prevprec { cout.precision(12) } ;
	cout << label << " = " << H << "\n" ;
	cout << std::setprecision (prevprec) ;
	}
    return H ;
    }

const Dvec& Numerics::eval_grad (bool print)		// Evaluate gradient vector
    {
    const auto&	use	{ eval_inuse (0, false) } ;
    const auto&	grad	{ global.data().grad } ;
    const auto&	Hterms	{ global.info().Hterms } ;
    const auto&	gens	{ global.info().gens[0] } ;
    int		neven	{ global.info().neven[0] } ;
    int		nterms	( Hterms.size() ) ;
    auto	recptr	{ grad.begin() } ;

    if (grad.entry().ncol != nterms ||
	grad.entry().nrow != neven)
	    {
	    cout << " ncol " << grad.entry().ncol << " nterms " << nterms << "\n" ;
	    cout << " nrow " << grad.entry().nrow << " neven " << neven << "\n" ;
	    gripe ("Need to (re)build gradient!") ;
	    }

    gradient.zeros(neven) ;
    for (int i(0) ; i < nterms ; ++i)
	{
	doub coeff { Hterms[i].coeff() } ;
	for (int j(0) ; j < neven ; ++j, ++recptr)
	    {
	    if (gens[j].active)
		{
		const Poly* ptr { &(*recptr) } ;
		const Poly* end { ptr++->end() } ;
		doub	    val ( 0.0 ) ;
		while (ptr < end) val += termvalue (ptr) ;
		gradient [j] += coeff * val ;
		}
	    }
	}
    if (recptr != grad.end()) fatal ("Inconsistent gradient record!") ;

    gradient = gradient (use) ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	gradient.raw_print ("Gradient =") ;
	cout << std::setprecision (prevprec) ;
	}
    return gradient ;
    }

const Dmtx& Numerics::eval_curv (uint repnum, int print)	// Evaluate T-even curvature
    {
    const auto&	use	{ eval_inuse (repnum, false) } ;
    const auto&	curv	{ global.data().curv[repnum] } ;
    const auto&	Hterms	{ global.info().Hterms } ;
    const auto&	gens	{ global.info().gens[repnum] } ;
    int		neven	{ global.info().neven[repnum] } ;
    int		ngen	( gens.size() ) ;
    int		nterms	( Hterms.size() ) ;
    auto	recptr	{ curv.begin() } ;

    if (curv.entry().nslice != nterms ||
	curv.entry().ncol != ngen  ||
	curv.entry().nrow != ngen) gripe ("Need to (re)build curvature!") ;

    curvature.zeros(neven, neven) ;
    for (int i(0) ; i < nterms ; ++i)
	{
	doub coeff { Hterms[i].coeff() } ;
	for (int j(0) ; j < ngen ; ++j)
	    {
	    for (int k(0) ; k < ngen ; ++k, ++recptr)
		{
		if (j < neven && gens[j].active &&
		    k < neven && gens[k].active)
		    {
		    const Poly*	ptr { &(*recptr) } ;
		    const Poly*	end { ptr++->end() } ;
		    doub	val ( 0.0 ) ;
		    while (ptr < end) val += termvalue (ptr) ;
		    curvature (j,k) += coeff * val ;
		    }
		}
	    }
	}
    if (recptr != curv.end()) fatal ("Inconsistent curvature record!") ;

    curvature = curvature (use,use) ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	if (print > 1)
	    {
	    cout << Rep::list[repnum].name ;
	    curvature.raw_print (" T-even curvature =") ;
	    }
	cout << Rep::list[repnum].name ;
	arma::sort (arma::eig_gen (curvature)).raw_print (" T-even curvature eigenvalues =") ;
	cout << std::setprecision (prevprec) ;
	}
    return curvature ;
    }

const Dmtx& Numerics::eval_metr (uint repnum, int print)	// Evaluate T-odd curvature
    {
    const auto&	use	{ eval_inuse (repnum, true) } ;
    const auto&	curv	{ global.data().curv[repnum] } ;
    const auto&	Hterms	{ global.info().Hterms } ;
    const auto&	gens	{ global.info().gens[repnum] } ;
    int		ngen	( global.info().gens[repnum].size() ) ;
    int		neven	{ global.info().neven[repnum] } ;
    int		nodd	{ ngen - neven } ;
    int		nterms	( Hterms.size() ) ;
    auto	recptr	{ curv.begin() } ;

    if (curv.entry().nslice != nterms ||
	curv.entry().ncol != ngen  ||
	curv.entry().nrow != ngen) gripe ("Need to (re)build curvature!") ;

    metric.zeros(nodd, nodd) ;
    for (int i(0) ; i < nterms ; ++i)
	{
	doub coeff { Hterms[i].coeff() } ;
	for (int j(0) ; j < ngen ; ++j)
	    {
	    for (int k(0) ; k < ngen ; ++k, ++recptr)
		{
		if (j < ngen && j >= neven && gens[j].active &&
		    k < ngen && k >= neven && gens[k].active)
		    {
		    const Poly* ptr { &(*recptr) } ;
		    const Poly* end { ptr++->end() } ;
		    doub	val ( 0.0 ) ;
		    while (ptr < end) val += termvalue (ptr) ;
		    metric (j-neven,k-neven) += coeff * val ;
		    }
		}
	    }
	}
    if (recptr != curv.end()) fatal ("Inconsistent curvature record!") ;

    metric = metric (use,use) ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	if (print > 1)
	    {
	    cout << Rep::list[repnum].name ;
	    metric.raw_print (" T-odd curvature =") ;
	    }
	cout << Rep::list[repnum].name ;
	arma::sort (arma::eig_gen (metric)).raw_print (" T-odd curvature eigenvalues =") ;
	cout << std::setprecision (prevprec) ;
	}
    return metric ;
    }

const Dmtx& Numerics::eval_lagr (uint repnum, int print)	// Evaluate Lagrange bracket matrix
    {
    const auto&	evens	{ eval_inuse (repnum, false) } ;
    const auto&	odds	{ eval_inuse (repnum, true)  } ;
    const auto&	lagr	{ global.data().lagr[repnum] } ;
    const auto&	gens	{ global.info().gens[repnum] } ;
    int		neven	{ global.info().neven[repnum] } ;
    int		ngen	( gens.size() ) ;
    int		nodd	{ ngen - neven } ;
    auto	recptr	{ lagr.begin() } ;

    if (lagr.entry().ncol != neven ||
	lagr.entry().nrow != nodd) gripe ("Need to build Lagrange matrix!") ;

    lagrange.zeros(neven, nodd) ;
    for (int i(0) ; i < neven ; ++i)
	{
	for (int j(neven) ; j < ngen ; ++j, ++recptr)
	    {
	    if (gens[i].active && gens[j].active)
		{
		const Poly* ptr { &(*recptr) } ;
		const Poly* end { ptr++->end() } ;
		doub	    val ( 0.0 ) ;
		while (ptr < end) val += termvalue (ptr) ;
		lagrange (i,j-neven) = val ;
		}
	    }
	}
    if (recptr != lagr.end()) fatal ("Inconsistent Lagrange bracket record!") ;

    lagrange = lagrange (evens,odds) ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	if (print > 1)
	    {
	    cout << Rep::list[repnum].name ;
	    lagrange.raw_print (" Lagrange bracket =") ;
	    }
	cout << Rep::list[repnum].name ;
	arma::svd(lagrange).raw_print (" Lagrange bracket singular values =") ;
	cout << std::setprecision (prevprec) ;
	}
    return lagrange ;
    }

doub Numerics::eval_delta (bool print)			// Evaluate delta vector
    {
    const auto&	use   { eval_inuse (0, false) } ;
    const auto&	grad  { eval_grad  () } ;
    const auto&	curv  { eval_curv  (0,0) } ;
    int		ngens { global.info().neven.front() } ;
    Dvec	del ;

    check_curv (curv) ;
    if (svdcut)	del = -arma::pinv (curv,svdcut) * grad ;
    else	del = -arma::solve (curv,grad) ;

    delta.zeros(ngens) ;
    delta(use) = del ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	delta.raw_print ("Delta =") ;
	cout << std::setprecision (prevprec) ;
	}
    return arma::norm (delta,"inf") ;
    }

const Cvec& Numerics::eval_spectra (string word, bool print)	// Evaluate particle spectrum
    {
    try { return eval_spectra (Rep::known (word), print) ; }
    catch (const exception& e) { gripe ("Unknown representation " + word) ; }
    }

const Cvec& Numerics::eval_spectra (uint repnum, bool print)	// Evaluate particle spectrum
    {
    auto&	metr	{ eval_metr (repnum) } ;
    auto&	lagr	{ eval_lagr (repnum) } ;
    Dmtx	curv	{ eval_curv (repnum) } ;
    Dmtx	inertia	{ lagr * metr.i() * lagr.t() } ;

    status.reset() ;
    if (global.symcurv) curv = (curv + curv.t()) / 2.0 ;
    check_curv (curv) ;

    if (!arma::eig_pair (spectrum, modes, curv, inertia))
	abort ("Cannot solve generalized eigensystem") ;
    if (curv_rpt (repnum)) cout << "\n" ;

    lastrep  = repnum ;
    spectrum = arma::sqrt (spectrum) ;
    if (!spectrum.has_nan()) spectrum = arma::sort (spectrum) ;
    if (print) Print::print_spectrum() ;
    return spectrum ;
    }

void Numerics::eval_geos (int printlim)			// Evaluate vev derivatives
    {
    eval_delta () ;
    Numerics::do_dvev (0.0, vev, dvev) ;

    if (printlim > dvev.size() || !printlim)
	printlim = dvev.size() ;
    cout << "d(Vev) = " ;
    for (int k(0) ; k < printlim ; ++k)
	{
	cout << (k % 5 ? ", " : "\n") << dvev[k] ;
	}
    cout << "\n" ;
    }

void Numerics::do_dvev (doub s, const Rvec& v, Rvec& dv)	// Evaluate vev derivs
    {
    const auto&	blab	{ Blab::level(Blab::NUMERICS) } ;
    const auto&	geos	{ global.data().geos } ;
    const auto&	bckt	{ global.info().bckt } ;
    long	offset	{ global.stage ? global.info(0).nobs : 0 } ;
    int		ngens	{ global.info().neven.front() } ;
    long	n	{ global.info().nobs } ;

    if (global.interrupt) return ;
    if (numerics.delta.n_elem != ngens) gripe ("Need to (re)evaluate delta!"); 
    if (geos[0].entry().nrow != ngens || !geos[0].entry().ncol)
	gripe ("Need to (re)build geodesic equations!") ;

    if (blab > 2) cout << format("do_dvev: s = {:.6f}, nvev = {}, ",s,v.n_elem) ;

    dv.zeros(n) ;
    numerics.dvev_buf = dv.memptr() ;

    if (v.memptr() != &numerics.vev[offset])
	{
	numerics.vev_tmp = numerics.vev ;
	if (global.stage)	numerics.vev_tmp.tail (v.n_elem) = v ;
	else			numerics.vev_tmp.head (v.n_elem) = v ;
	numerics.vev_buf = numerics.vev_tmp.memptr() ;
	}
    else numerics.vev_buf = numerics.vev.memptr() ;

    TASK_ARENA (global.maxthread, bckt,
	FOR_EACH (bckt.begin(), bckt.end(), do_dvev_bckt)) ;

    if (blab > 2)
	{
	doub	maxv  (0.0) ;
	doub	maxdv (0.0) ;
	numb	maxi  (0) ;
	numb	maxdi (0) ;
	auto	nvev { v.n_elem } ;

	for (int indx(global.stage ? 0 : 1) ; indx < nvev ; ++indx)
	    {
	    if (std::abs(v[indx]) > maxv)
		{
		maxv = std::abs(v[indx]) ;
		maxi = indx ;
		}
	    if (std::abs(dv[indx]) > maxdv)
		{
		maxdv = std::abs(dv[indx]) ;
		maxdi = indx ;
		}
	    }
	cout << format("maxv #{}={:.5g}, ", maxi,  maxv)  ;
	cout << format("maxdv #{}={:.5g}\n", maxdi, maxdv) << flush ;
	}
    }

void Numerics::do_dvev_bckt (const numb3& bucket)		// Evaluate dvev bucket
    {
    if (global.interrupt) return ;

    numb	bcktnum	{ bucket[0] } ;
    numb	first	{ bucket[1] } ;
    numb	last	{ bucket[2] } ;
    long	offset	{ global.stage ? global.info(0).nobs : 0 } ;
    const auto&	delta	{ numerics.delta } ;
    const auto&	geos	{ global.data().geos[bcktnum] } ;
    real*	dv 	{ numerics.dvev_buf } ;
    const real*	v	{ numerics.vev_buf } ;
    int		ngens	( delta.size() ) ;

    if (geos.entry().ncol != last - first + 1 ||
	geos.entry().nrow != ngens)
	{
	cout << "geo bckt " << bcktnum << " first " << first << " last " << last
	     << " ncol " << geos.entry().ncol << " expected " << last - first + 1
	     << " gens nrow " << geos.entry().nrow << " expected " << ngens << "\n" << flush ;
	gripe ("\nNeed to (re)build geodesic equations!") ;
	}

    if (global.geoswap) Save::read_geo_bckt (bcktnum) ;
    auto recptr	{ geos.begin() } ;
    for (int i(first) ; i <= last ; ++i)
	{
	for (int j(0) ; j < ngens ; ++j, ++recptr)
	    {
	    if (global.interrupt) return ;
	    if (doub coef { delta[j] })
		{
		const Poly* ptr  { &(*recptr) } ;
		const Poly* end  { ptr++->end() } ;
		doub	    val  ( 0.0 ) ;
		while (ptr < end) val += termvalue (ptr,v) ;
		dv [i - offset] += coef * val ;
		}
	    }
	}
    if (recptr != geos.end())
	fatal ("Inconsistent geodesic record") ;
    if (global.geoswap) Save::read_geo_bckt (-bcktnum-1) ;
    }

doub Numerics::err_norm (const Rvec& err, const Rvec& y)	// ODE error vector norm
    //
    //	Evaluate ODE error norm
    //
    {
    if (global.interrupt) return 0 ;

    const auto&	blab	{ Blab::level(Blab::NUMERICS) } ;
    doub	norm	{ arma::norm (err,"inf") } ;
//    auto&	eps	{ numerics.odetol } ;
//    doub	norm	{ arma::norm (err / (arma::abs (y) + eps),"inf") } ;
//    doub	norm	{ arma::norm (err,2) / sqrt (err.n_elem) } ;
    if (blab > 2) cout << "err_norm: " << norm << "\n" << flush ;
    return norm ;
    }

void Numerics::status_rpt (uint iters, uint steps)
    {
    auto&	okneg	{ global.oknegeig } ;
    bool	euclidF { global.stage && theory.euclid } ;
    doub&	mostneg { status.negeig[0] } ;
    auto&	maxi	{ status.maxloopi  } ;
    auto&	maxv	{ status.maxloopv  } ;
    bool	unphys  { !global.stage && check_loops() } ;

    if (!curv_rpt() && !unphys) cout << " ok" ;
    else if (unphys) cout << format(" unphys (#{}={:3.1f})", maxi, maxv) ;
    cout << " " << iters << "/" << steps << " iters/steps\n" << flush ;

    if (!euclidF && !okneg && mostneg < -svdcut)
	abort ("Negative curvature eigenvalues") ;
    }

bool Numerics::open_MMA ()				// Open MMA output file
    {
    const auto&	MMAdir		{ global.MMAdir  } ;
    auto&	MMApath		{ global.info().MMAfile.path } ;
    auto&	MMAstream	{ global.info().MMAfile.stream } ;
    string	file		{ global.mk_filename("m") } ;
    string	path		{ global.addsubdir (MMAdir) + file } ;
    auto	mode		{ std::ios::out } ;
    string	okmsg		{ "Writing results to " } ;
    ofstream	stream ;

    if (path != MMApath) MMAstream.close() ;
    if (!MMAstream.is_open())
	{
	if (global.info().MMAfile.append)
	    {
	    mode |= std::ios::app ;
	    okmsg = "Appending results to " ;
	    }
	stream.open (path, mode) ;
	if (stream.is_open())
	    {
	    cout << okmsg << path << "\n" ;
	    MMApath	  = std::move (path) ;
	    MMAstream = std::move (stream) ;

	    string sep {""} ;
	    MMAstream << "coupling = { " ;
	    for (const auto& c : Coupling::list)
		{
		if (c.stage > global.stage) continue ;
		MMAstream << sep << '"' << c << '"' ;
		sep = ", " ;
		}
	    MMAstream << " } ;\n" ;
	    }
	else gripe ("Cannot write MMA results file " + path) ;
	}
    return MMAstream.good() ;
    }

void Numerics::write_data ()				// Write data to MMAfile
    {
    if (!open_MMA()) gripe ("Cannot write MMA results file!") ;

    string	HorF	{ theory.euclid ? "F" : "H" } ;
    auto&	stream	{ global.info().MMAfile.stream } ;

    data_write (stream, HorF, eval_ham()) ;

    for (auto& [i,o] : global.info().MMAfile.obs)
	{
	data_write (stream, o.print(), vev[i]) ;
	}

    if (!theory.euclid)
	{
	for (int i(0) ; i < Rep::list.size() ; ++i)
	    {
	    if (!built_rep(i)) continue ;
	    data_write (stream, Rep::list[i].name, eval_spectra(i)) ;
	    }
	}
    stream << flush ;
    }

bool Numerics::check_loops ()					// Loop vevs < 1?
    {
    if (ObsList::swapped) return false ;
    
    const auto&	obslist { ObsList::obs } ;
    auto 	beg	{ obslist.begin() } ;
    auto 	end	{ obslist.end()   } ;
    doub	maxv	(0.0) ;
    uint	maxi	(0) ;

    for (auto ptr { beg } ; ptr < end ; ++ptr)
	{
	int	indx ( ptr - beg ) ;
	auto&	obs  { obslist(indx) } ;
	if (obs.is_Loop() && std::abs(vev[indx]) > maxv)
	    {
	    maxv = std::abs(vev[indx]) ;
	    maxi = indx ;
	    }
	}
    status.maxloopv = vev[maxi] ;
    status.maxloopi = maxi ;
    return maxv > 1.0 ;
    }

bool Numerics::check_curv (const Dmtx& curv)			// OK curvature eigs?
    {
    if (curv.n_rows) 
	{
	Cvec	eigs  { arma::eig_gen(curv) } ;
	Uvec	reals { arma::sort_index (arma::real (eigs),"ascent")  } ;
	Uvec	imags { arma::sort_index (arma::imag (eigs),"descend") } ;
	bool	negok { theory.euclid && global.stage } ;
	uint	minr  ( reals[0] ) ;
	uint	maxi  ( imags[0] ) ;
	int	ncmplx (0) ;
	int	nneg   (0) ;

	if (eigs[minr].real() < 0)
	    {
	    status.warn = true ;
	    for (auto i : reals)
		{
		if (eigs[i].real() >= 0) break ;
		if (nneg == status.negeig.size()) break ;
		status.negeig[nneg++] = eigs[i].real() ;
		}
	    }
	if (eigs[maxi].imag() > 0)
	    {
	    status.warn = true ;
	    for (auto i : imags)
		{
		if (eigs[i].imag() <= 0) break ;
		if (ncmplx == status.cmplxeig.size()) break ;
		status.cmplxeig[ncmplx++] = eigs[i] ;
		}
	    }
	}
    return status.warn ;
    }

bool Numerics::curv_rpt (int repnum)
    {
    if (status.warn)
	{
	if (repnum >= 0) cout << Rep::list[repnum].name << ": " ;
	if (status.negeig[0] < 0)
	    {
	    char prefix {'('} ;
	    cout << " negeigs " ;
	    for (auto& z : status.negeig)
		{
		if (!z) break ;
		cout << format("{}{:4.2f}", prefix, z) ;
		prefix = ',' ;
		}
	    cout << ')' ;
	    }
	if (status.cmplxeig[0].imag())
	    {
	    char prefix {'('} ;
	    cout << " cmplxeigs " ;
	    for (auto& z : status.cmplxeig)
		{
		if (!z.imag()) break ;
		cout << format("{}{:3.1f}{:+4.2f}i", prefix, z.real(), z.imag()) ;
		prefix = ',' ;
		}
	    cout << ')' ;
	    }
	status.reset() ;
	return true ;
	}
    return false ;
    }

bool Numerics::built_rep (uint repnum)				// Is irrep built?
    {
    const auto&	curv	{ global.data().curv[repnum] } ;
    const auto&	lagr	{ global.data().lagr[repnum] } ;
    int		Hterms	( global.info().Hterms.size() ) ;
    auto	ngens	{ global.info().gens[repnum].size() } ;
    auto	neven	{ global.info().neven[repnum] } ;
    auto	nodd	{ ngens - neven } ;

    return  lagr.entry().nrow == neven && lagr.entry().ncol == nodd  &&
	    curv.entry().nrow == ngens && curv.entry().ncol == ngens && 
	    curv.entry().nslice == Hterms && ngens > 0 ;
    }

void Numerics::data_write (ofstream& out, const string s, doub v)	// Save value
    {
    string c { Coupling::values() } ;
    out << s << "[" << c << "] = " << MMAform(v) << " ;\n" ;
    }

void Numerics::data_write (ofstream& out, const string s, const Rvec& vec) // Save vector
    {
    string c { Coupling::values() } ;
    if (vec.size())
	{
	string delim { "{ " } ;
	out << s << "[" << c << "] = " ;
	for (const auto& v : vec)
	    {
	    out << delim << MMAform(v) ;
	    delim = ", " ;
	    }
	out << " } ;\n" ;
	}
    }

void Numerics::data_write (ofstream& out, const string s, const Cvec& vec) // Save vector
    {
    string c { Coupling::values() } ;
    if (vec.size())
	{
	string delim { "{ " } ;
	out << s << "[" << c << "] = " ;
	for (const auto& v : vec)
	    {
	    out << delim << MMAform(v.real()) ;
	    if (v.imag()) out << " + " << MMAform(v.imag()) << " I" ;
	    delim = ", " ;
	    }
	out << " } ;\n" ;
	}
    }

string Numerics::MMAform (doub x)			// Convert to MMA input form
    {
    std::stringstream buf ;
    buf << std::setprecision(12) << x ;
    string	s   { buf.str() } ;
    auto	pos { s.find ("e") } ;
    if (pos != s.npos) { s.replace(pos, 1, "*^") ; }
    return s ;
    }

void Numerics::initialize (int stage)			// Initialize expectation values
    {
    const auto&	blab	{ Blab::level(Blab::NUMERICS) } ;
    const auto& nobsG	{ global.info(0).nobs } ;
    const auto& nobsF	{ global.info(1).nobs } ;

    vev.resize(nobsG + nobsF) ;
    vev[0] = 1.0 ;

    if (stage == 0)				// gauge vev's
	{
	vev.zeros() ;
	vev[0] = 1.0 ;
	if (blab) cout << "(Re)initialized gauge vev's\n" ;
	}
    if (nobsF)					// fermion vev's
	{
	vev.tail(nobsF).zeros() ;

	char8	mass { "mass" } ;
	int	indx { Coupling::indx (mass) } ;
	auto&	m    { Coupling::list[indx].value } ;
	doub	condensate ;

	m = Coupling::dfltval ;
	if (theory.euclid)	condensate = -1.0 / m ;
	else			condensate = m > 0 ? -0.5 : 0.5 ;

	for (const auto& [indx_f,indx_g] : ObsList::fermiinit)
	    {
	    vev[indx_f] = vev[indx_g] * condensate ;
	    }
	if (blab) cout << "(Re)initialized fermion vev's\n" ;
	}
    global.fermivev = false ;
    }
