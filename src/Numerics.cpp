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

void Numerics::do_minimize ()				// Do minimization
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

int Numerics::do_step (doub tol)			// Do geodesic integration step
    {
    if (global.interrupt) return 0 ;

    uint	blab	{ Blab::level(Blab::NUMERICS) } ;
    Ode		ode	{ do_dvev, err_norm, odetol, rk, maxode } ;
    doub	dnorm	{ eval_delta() } ;
    uint	n	{ global.info().nobs } ;
    doub	s (0) ;
    bool	ok ;

    if (blab > 1) cout << "  |delta| = " << dnorm << ", " ;
    if (blab > 2) raw_print (delta, "\n   delta =") ;

    if (tol && dnorm <= tol) return 0 ;			// converged?

    if (n != vev.size())				// integrate vev subvec
	{
	real*	vevptr { &vev [global.stage ? global.info(0).nobs : 0] } ; 
	Rvec	subvev { aliasvec (vevptr, n) } ;
	ok = ode.integrate (s, 1.0, subvev) ;
	if (memptr (subvev) != vevptr)			// just for Eigen (sigh)
	    {
	    if (global.stage)	vev.tail (n) = subvev ;
	    else		vev.head (n) = subvev ;
	    }
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

    set_size (inuse, end-beg) ;
    for (int i(beg) ; i < end ; ++i)
	if (gens[i].active) inuse (n++) = i - beg ;
    resize (inuse,n) ;
    return inuse ;
    }

doub Numerics::eval_ham (bool print)			// Evaluate Hamiltonian/free energy
    {
    if (!global.maxord()) gripe ("Make observables first!") ;

    H = 0 ;
    for (const auto& Hterm : global.info().Hterms)
	{
	doub val (0.0) ;
	for (const auto& term : Hterm.cpoly)
	    {
	    val += termvalue (term) ;
	    }
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
    int		Hterm   { -1 } ;

    if (grad.entry().ncol != nterms ||
	grad.entry().nrow != neven)
	    {
	    cout << " ncol " << grad.entry().ncol << " nterms " << nterms << "\n" ;
	    cout << " nrow " << grad.entry().nrow << " neven " << neven << "\n" ;
	    gripe ("Need to (re)build gradient!") ;
	    }

    set_zero (gradient, neven) ;
    for (const auto& poly : grad)
	{
	const GradHdr& info ( poly ) ;
	if (gens[info.gen].active)
	    {
	    doub val ( 0.0 ) ;
	    for (const auto& term : poly) val += termvalue (term) ;
	    if (Hterm != info.Hterm)
		Hterm  = info.Hterm ;
	    gradient [info.gen] += Hterms[Hterm].coeff() * val ;
	    }
	}
    gradient = gradient (use) ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	raw_print (gradient,"Gradient =") ;
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
    int		Hterm   { -1 } ;

    if (curv.entry().nslice != nterms ||
	curv.entry().ncol != ngen  ||
	curv.entry().nrow != ngen) gripe ("Need to (re)build curvature!") ;

    set_zero (curvature, neven, neven) ;
    for (const auto& poly : curv)
	{
	const CurvHdr& info ( poly ) ;
	if (info.gen1 < neven && gens[info.gen1].active &&
	    info.gen2 < neven && gens[info.gen2].active)
	    {
	    doub val ( 0.0 ) ;
	    for (const auto& term : poly) val += termvalue (term) ;
	    if (Hterm != info.Hterm)
		Hterm  = info.Hterm ;
	    curvature (info.gen1,info.gen2) += Hterms[Hterm].coeff() * val ;
	    }
	}
    curvature = curvature (use,use) ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	if (print > 1)
	    {
	    cout << Rep::list[repnum].name ;
	    raw_print (curvature," T-even curvature =") ;
	    }
	cout << Rep::list[repnum].name ;
	raw_print (sort (eig_gen (curvature))," T-even curvature eigenvalues =") ;
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
    int		Hterm   { -1 } ;
    doub	coeff ;

    if (curv.entry().nslice != nterms ||
	curv.entry().ncol != ngen  ||
	curv.entry().nrow != ngen) gripe ("Need to (re)build curvature!") ;

    set_zero (metric, nodd, nodd) ;
    for (const auto& poly : curv)
	{
	const CurvHdr&	info ( poly ) ;
	if (info.gen1 < ngen && info.gen1 >= neven && gens[info.gen1].active &&
	    info.gen2 < ngen && info.gen2 >= neven && gens[info.gen2].active)
	    {
	    doub val  ( 0.0 ) ;
	    for (const auto& term : poly) val += termvalue (term) ;
	    if (Hterm != info.Hterm)
		Hterm  = info.Hterm ;
	    metric (info.gen1-neven,info.gen2-neven) += Hterms[Hterm].coeff() * val ;
	    }
	}
    metric = metric (use,use) ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	if (print > 1)
	    {
	    cout << Rep::list[repnum].name ;
	    raw_print (metric," T-odd curvature =") ;
	    }
	cout << Rep::list[repnum].name ;
	raw_print (sort (eig_gen (metric))," T-odd curvature eigenvalues =") ;
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

    if (lagr.entry().ncol != ngen ||
	lagr.entry().nrow != ngen) gripe ("Need to build Lagrange matrix!") ;

    set_zero (lagrange, neven, nodd) ;
    for (const auto& poly : lagr)
	{
	const LagrHdr& info ( poly ) ;
	if (info.gen1 < neven && gens[info.gen1].active &&
	    info.gen2 < ngen  && gens[info.gen2].active && info.gen2 >= neven)
	    {
	    doub val ( 0.0 ) ;
	    for (const auto& term : poly) val += termvalue (term) ;
	    lagrange (info.gen1,info.gen2-neven) = val ;
	    }
	}
    lagrange = lagrange (evens,odds) ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	if (print > 1)
	    {
	    cout << Rep::list[repnum].name ;
	    raw_print (lagrange," Lagrange bracket =") ;
	    }
	cout << Rep::list[repnum].name ;
	raw_print (svd(lagrange)," Lagrange bracket singular values =") ;
	cout << std::setprecision (prevprec) ;
	}
    return lagrange ;
    }

doub Numerics::eval_delta (bool print)		// Evaluate delta vector
    {
    const auto&	use   { eval_inuse (0, false) } ;
    const auto&	grad  { eval_grad  () } ;
    const auto&	curv  { eval_curv  (0,0) } ;
    int		ngens { global.info().neven.front() } ;
    Dvec	del ;

    check_curv (curv) ;
    if (svdcut)	del = -pinv (curv,svdcut) * grad ;
    else	del = -linsolve (curv,grad) ;

    set_zero (delta, ngens) ;
    delta(use) = del ;
    if (print)
	{
	auto prevprec { cout.precision(12) } ;
	raw_print (delta, "Delta =") ;
	cout << std::setprecision (prevprec) ;
	}
    return infnorm (delta) ;
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
    Dmtx	inertia	{ lagr * inv(metr) * transpose(lagr) } ;

    status.reset() ;
    if (global.symcurv) curv = (curv + transpose(curv)) / 2.0 ;
    check_curv (curv) ;

    if (!eig_pair (spectrum, modes, curv, inertia))
	abort ("Cannot solve generalized eigensystem") ;
    if (curv_rpt (repnum)) cout << "\n" ;

    lastrep  = repnum ;
    spectrum = sqrt (spectrum) ;
    if (!has_nan (spectrum)) spectrum = sort (spectrum) ;
    if (print) Print::print_spectrum() ;
    return spectrum ;
    }

void Numerics::eval_geos (int printlim)			// Evaluate observable derivatives
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
    uint	blab	{ Blab::level(Blab::NUMERICS) } ;
    uint	offset	{ global.stage ? global.info(0).nobs : 0 } ;
    const auto&	geos	{ global.data().geos } ;
    const auto&	bckt	{ global.info().bckt } ;
    int		ngens	{ global.info().neven.front() } ;
    uint	n	{ global.info().nobs } ;

//    std::sort (bckt.begin(), bckt.end(), [geos](const uint3& a, const uint3& b)
//	{ return geos[a[0]].entry().filepos < geos[b[0]].entry().filepos ; }) ;

    if (global.interrupt) return ;
    if (nelem(numerics.delta) != ngens) gripe ("Need to (re)evaluate delta!"); 
    if (geos[0].entry().nrow != ngens || !geos[0].entry().ncol)
	gripe ("Need to (re)build geodesic equations!") ;

    if (blab > 2) cout << format("do_dvev: s = {:.6f}, nvev = {}, ",s,nelem(v)) ;

    set_zero (dv, n) ;
    numerics.dvev_buf = memptr(dv) ;

    if (memptr(v) != &numerics.vev[offset])
	{
	numerics.vev_tmp = numerics.vev ;
	if (global.stage)	numerics.vev_tmp.tail (nelem(v)) = v ;
	else			numerics.vev_tmp.head (nelem(v)) = v ;
	numerics.vev_buf = memptr(numerics.vev_tmp) ;
	}
    else numerics.vev_buf = memptr(numerics.vev) ;

    TASK_ARENA (global.maxthread, bckt,
	FOR_EACH (bckt.begin(), bckt.end(), do_dvev_bckt)) ;

    if (blab > 2)
	{
	doub	maxv  (0.0) ;
	doub	maxdv (0.0) ;
	uint	maxi  (0) ;
	uint	maxdi (0) ;
	auto	nvev { nelem(v) } ;

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

void Numerics::do_dvev_bckt (const uint3& bucket)		// Evaluate dvev bucket
    {
    if (global.interrupt) return ;

    uint	bcktnum	{ bucket[0] } ;
    uint	first	{ bucket[1] } ;
    uint	last	{ bucket[2] } ;
    uint	offset	{ global.stage ? global.info(0).nobs : 0 } ;
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
    for (const auto& poly : geos)
	{
	if (global.interrupt) return ;

	const GeoHdr&	info ( poly ) ;
	real		val  ( 0.0 ) ;
	doub		coef { delta[info.gen] } ;

	if (info.indx < first || info.indx > last)
	    fatal (format ("Inconsistent geo record: bucket {} indx {} first {} last {}!",
		    bcktnum, info.indx, first, last)) ; 

	if (global.interrupt) break ;
	if (coef)
	    {
	    for (const auto& term : poly) val += termvalue (term,v) ;
	    dv [info.indx - offset] += coef * val ;
	    }
	}
    if (global.geoswap) Save::read_geo_bckt (-bcktnum-1) ;
    }

doub Numerics::err_norm (const Rvec& err, const Rvec& y)	// ODE error vector norm
    //
    //	Evaluate ODE error norm
    //
    {
    if (global.interrupt) return 0 ;

    uint	blab { Blab::level(Blab::NUMERICS) } ;
//    auto&	eps  { numerics.odetol } ;
//    ulong	maxi { index_max (abs (err)) } ;
//    doub	norm { l2norm  (err) / sqrt (nelem(err)) } ;
//    doub	norm { infnorm (err / (abs (err) + eps)) } ;
    doub	norm { infnorm (err) } ;
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
    auto&	MMApath		{ global.info().MMApath } ;
    auto&	MMAstream	{ global.info().MMAstream } ;
    string	file		{ global.mk_filename("m") } ;
    string	path		{ global.addsubdir (MMAdir) + file } ;
    auto	mode		{ std::ios::out } ;
    string	okmsg		{ "Writing results to " } ;
    ofstream	stream ;

    if (path != MMApath) MMAstream.close() ;
    if (!MMAstream.is_open())
	{
	if (global.MMAappend)
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
    auto&	stream	{ global.info().MMAstream } ;

    data_write (stream, HorF, eval_ham()) ;

    for (auto& [i,o] : global.info().MMAobs)
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

bool Numerics::check_loops ()						// Loop vevs < 1?
    {
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

bool Numerics::check_curv (const Dmtx& curv)				// OK curvature eigs?
    {
    if (nrow(curv)) 
	{
	Cvec	eigs  { eig_gen(curv) } ;
	Uvec	reals { sort_index (realpart (eigs)) } ;
	Uvec	imags { sort_index (imagpart (eigs), false) } ;
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

bool Numerics::built_rep (uint repnum)					// Is irrep built?
    {
    const auto&	curv	{ global.data().curv[repnum] } ;
    const auto&	lagr	{ global.data().lagr[repnum] } ;
    int		ngens	( global.info().gens[repnum].size() ) ;
    int		Hterms	( global.info().Hterms.size() ) ;

    return  lagr.entry().nrow == ngens && lagr.entry().ncol == ngens &&
	    curv.entry().nrow == ngens && curv.entry().ncol == ngens && 
	    curv.entry().nslice == Hterms && ngens > 0 ;
    }

void Numerics::data_write (ofstream& out, const string s, doub v)		// Save data value
    {
    string c { Coupling::values() } ;
    out << s << "[" << c << "] = " << MMAform(v) << " ;\n" ;
    }

void Numerics::data_write (ofstream& out, const string s, const Rvec& vec)	// Save data vector
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

void Numerics::data_write (ofstream& out, const string s, const Cvec& vec)	// Save data vector
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

string Numerics::MMAform (doub x)				// Convert to MMA input form
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
    uint	blab	{ Blab::level(Blab::NUMERICS) } ;
    const auto& nobsG	{ global.info(0).nobs } ;
    const auto& nobsF	{ global.info(1).nobs } ;
    resize (vev, nobsG + nobsF) ;
    vev[0] = 1.0 ;

    if (stage == 0)				// gauge vev's
	{
	set_zero (vev) ;
	vev[0] = 1.0 ;
	if (blab) cout << "(Re)initialized gauge vev's\n" ;
	}
    if (nobsF)					// fermion vev's
	{
	set_zero (vev.tail (nobsF)) ;

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
