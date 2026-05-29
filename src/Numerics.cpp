#include "Numerics.h"
#include "Global.h"
#include "Rep.h"
#include "Data.h"
#include "Save.h"
#include "Blab.h"
#include "Print.h"
#include "Gripe.h"
#include <cmath>

void Numerics::do_flow (int indx, doub v0, doub v1, doub inc)	// Flow coupling
    {
    auto&	coupling { Coupling::list[indx] } ;
    doub&	value    { coupling.value } ;
    auto	prevprec { cout.precision(12) } ;

    if (!inc) gripe ("Must have non-zero coupling increment!") ;
    if ((v1 - v0)/inc < 0) inc *= -1 ;

    cout << std::scientific ;
    for (value = v0 ;;)
	{
	cout << coupling << format(" = {:6.3f}:", value) << flush ;
	try {
	    do_minimize () ;
	    write_data (value) ;
	    }
	catch (const Abort& e)
	    {
	    cout.precision (prevprec) ;
	    cout << std::defaultfloat << e.what() << "\n" ;
	    abort (format ("Aborting flow at {} = {:.6f}", coupling.data(), value)) ;
	    }
	if (value == v1)				break ;
	else if (abs(value - v1) > abs(1.6 * inc))	value += inc ;
	else if (abs(value - v1) > abs(inc))		value  = (value + v1)/2 ;
	else						value  = v1 ;
	}
    cout << std::setprecision (prevprec) << std::defaultfloat ;
    }

void Numerics::do_minimize ()				// Do minimization
    {
    int  delta_s ;
    doub tol   (0) ;
    uint iters (0) ;
    uint steps (0) ;
    ++stats.tries ;
    status.reset() ;

    do  {
	if (global.interrupt || ++iters > minmax) break ;
	delta_s = do_step (tol) ;
	steps += abs(delta_s) ;
	tol = mintol ;
	} while (delta_s > 0) ;

    status_rpt (iters,steps) ;
    if (iters > minmax)
	{
	++stats.fails ;
	abort ("Minimization failed, max iterations exceeded") ;
	}
    else if (delta_s < 0)
	{
	++stats.fails ;
	abort ("ODE integration failed!") ;
	}
    if (global.autosave && !global.interrupt) Save::save_vev () ;
    }

int Numerics::do_step (doub tol)			// Do geodesic integration step
    {
    if (global.interrupt) return 0 ;

    uint	nvev  { global.info().nobs } ;
    uint	blab  { Blab::blablevel[BLAB::NUMERICS] } ;
    Ode		ode   { do_dvev, err_norm, odetol, rk, odemax } ;
    doub	dnorm { eval_delta() } ;
    doub	s (0) ;
    bool	ok ;

    if (blab > 1) cout << "  |delta| = " << dnorm << ", " ;
    if (blab > 2) cout << "\n   delta = \n" << delta << flush ;

    if (tol && dnorm <= tol) return 0 ;			// converged?

    if (nvev != global.nobs())				// integrate vev subvec
	{
	real*	vevptr { &vev [global.stage ? global.nobsG() : 0] } ; 
	Rvec	subvev { aliasvec (vevptr, nvev) } ;
	ok = ode.integrate (s, 1.0, subvev) ;
	}
    else ok = ode.integrate (s, 1.0, vev) ;

    if (blab > 1) cout << ode.steps << " step(s) " << ode.rejects << " rejects "
		       << (theory.euclid ? "F = " : "H = ") << eval_H() << "\n" << flush ;

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
	{
	if (gens[i].active) inuse (n++) = i - beg ;
	}
    inuse.resize(n) ;
    return inuse ;
    }

const Dvec& Numerics::eval_grad (bool print)		// Evaluate gradient vector
    {
    const auto& coup	{ Coupling::list } ;
    const auto&	use	{ eval_inuse (0, false) } ;
    const auto& grad	{ global.data().grad } ;
    const auto& Hterms	{ global.info().Hterms } ;
    const auto& gens	{ global.info().gens[0] } ;
    int		ngens	{ global.info().neven[0] } ;
    int		nterms	( Hterms.size() ) ;
    int		Hterm   { -1 } ;

    if (grad.entry().ncol != nterms ||
	grad.entry().nrow != ngens) gripe ("Need to (re)build gradient!") ;

    set_zero (gradient, ngens) ;
    for (const auto& poly : grad)
	{
	const GradHdr&	info ( poly ) ;
	if (gens[info.gen].active)
	    {
	    doub		val  ( 0.0 ) ;
	    for (const auto& term : poly) val += termvalue (term) ;
	    if (Hterm != info.Hterm)
		Hterm  = info.Hterm ;
	    gradient [info.gen] += Hterms[Hterm].coeff() * val ;
	    }
	}
    gradient = gradient (use) ;
    if (print) cout << "Gradient = \n" << gradient ;
    return gradient ;
    }

const Dmtx& Numerics::eval_curv (int repnum, int print)	// Evaluate T-even curvature
    {
    const auto& coup	{ Coupling::list } ;
    const auto&	use	{ eval_inuse (repnum, false) } ;
    const auto& curv	{ global.data().curv[repnum] } ;
    const auto& Hterms	{ global.info().Hterms } ;
    const auto& gens	{ global.info().gens[repnum] } ;
    int		ngen	( global.info().gens[repnum].size() ) ;
    int		neven	{ global.info().neven[repnum] } ;
    int		nterms	( Hterms.size() ) ;
    int		Hterm   { -1 } ;
    doub	coeff ;

    if (curv.entry().nslice != nterms ||
	curv.entry().ncol != ngen  ||
	curv.entry().nrow != ngen) gripe ("Need to (re)build curvature!") ;

    set_zero (curvature, neven, neven) ;
    for (const auto& poly : curv)
	{
	const CurvHdr&	info ( poly ) ;
	if (info.gen1 < neven && gens[info.gen1].active &&
	    info.gen2 < neven && gens[info.gen2].active)
	    {
	    doub val  ( 0.0 ) ;
	    for (const auto& term : poly) val += termvalue (term) ;
	    if (Hterm != info.Hterm)
		Hterm  = info.Hterm ;
	    curvature (info.gen1,info.gen2) += Hterms[Hterm].coeff() * val ;
	    }
	}
    curvature = curvature (use,use) ;
    if (print > 1)
	{
	cout << Rep::list[repnum].name << " T-even curvature = \n" ;
	raw_print (curvature, cout) ;
	}
    if (print)
	{
	cout << Rep::list[repnum].name << " T-even curvature eigenvalues = \n"
	     << sort (eig_gen (curvature)) ;
	}
    return curvature ;
    }

const Dmtx& Numerics::eval_metr (int repnum, int print)	// Evaluate T-odd curvature
    {
    const auto& coup	{ Coupling::list } ;
    const auto&	use	{ eval_inuse (repnum, true) } ;
    const auto& curv	{ global.data().curv[repnum] } ;
    const auto& Hterms	{ global.info().Hterms } ;
    const auto& gens	{ global.info().gens[repnum] } ;
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
    if (print > 1)
	{
	cout << Rep::list[repnum].name << " T-odd curvature = \n" << metric ;
	}
    if (print)
	{
	cout << Rep::list[repnum].name << " T-odd curvature eigenvalues = \n"
	     << sort (eig_gen (metric)) ;
	}
    return metric ;
    }

const Dmtx& Numerics::eval_lagr (uint repnum, bool print)	// Evaluate Lagrange bracket matrix
    {
    const auto&	evens	{ eval_inuse (repnum, false) } ;
    const auto&	odds	{ eval_inuse (repnum, true)  } ;
    const auto& lagr	{ global.data().lagr[repnum] } ;
    const auto& gens	{ global.info().gens[repnum] } ;
    int		ngen	( global.info().gens[repnum].size() ) ;
    int		neven	{ global.info().neven[repnum] } ;
    int		nodd	{ ngen - neven } ;

    if (lagr.entry().ncol != ngen ||
	lagr.entry().nrow != ngen) gripe ("Need to build Lagrange matrix!") ;

    set_zero (lagrange, neven, nodd) ;
    for (const auto& poly : lagr)
	{
	const LagrHdr&	info ( poly ) ;
	if (info.gen1 < neven && gens[info.gen1].active &&
	    info.gen2 < ngen  && gens[info.gen2].active && info.gen2 >= neven)
	    {
	    doub	val  ( 0.0 ) ;
	    for (const auto& term : poly) val += termvalue (term) ;
	    lagrange (info.gen1,info.gen2-neven) = val ;
	    }
	}
    lagrange = lagrange(evens,odds) ;
    if (print)
	{
	cout << Rep::list[repnum].name << " Lagrange bracket = \n" << lagrange ;
	cout << Rep::list[repnum].name << " Lagrange bracket singular values = \n"
	     << svd (lagrange) ;
	}
    return lagrange ;
    }

doub Numerics::eval_H (bool print)			// Evaluate Hamiltonian/free energy
    {
    auto	label	{ theory.euclid ? "Free energy" : "Hamiltonian" } ;
    const auto& Hterms	{ global.info().Hterms } ;
    const auto& coup	{ Coupling::list } ;

    H = 0 ;
    for (const auto& Hterm : Hterms)
	{
	doub	val (0.0) ;
	for (const auto& term : Hterm.cpoly) val += termvalue (term) ;
	H += Hterm.coeff() * val ;
	}
    if (print)
	{
	auto	prev { cout.precision(12) } ;
	cout << label << " = " << H << std::setprecision (prev) << "\n" ;
	}
    return H ;
    }

doub Numerics::eval_delta (bool print)		// Evaluate delta vector
    {
    const auto&	use   { eval_inuse (0, false) } ;
    const auto&	grad  { eval_grad () } ;
    const auto&	curv  { eval_curv (0,0) } ;
    int		ngens { global.info().neven.front() } ;
    Dvec	del ;

    check_curv (curv) ;
    if (svdcut)	del = -pinv (curv,svdcut) * grad ;
    else	del = -linsolve (curv, grad) ;

    set_zero (delta, ngens) ;
    delta(use) = del ;
    if (print) cout << "Delta = \n" << delta ;
    return infnorm (delta) ;
    }

const Cvec& Numerics::eval_spectra (string word, bool print) // Evaluate particle spectrum
    {
    try { return eval_spectra (Rep::known (word), print) ; }
    catch (const exception& e) { gripe ("Unknown representation " + word) ; }
    }

const Cvec& Numerics::eval_spectra (int repnum, bool print)	// Evaluate particle spectrum
    {
    auto&	metr	{ eval_metr  (repnum) } ;
    auto&	lagr	{ eval_lagr  (repnum) } ;
    Dmtx	curv	{ eval_curv  (repnum) } ;
    Dmtx	inertia { lagr * inv(metr) * transpose(lagr) } ;

    status.reset() ;
    if (global.symcurv) curv = (curv + transpose(curv)) / 2.0 ;
    check_curv (curv) ;

    if (!eig_pair (spectrum, modes, curv, inertia))
	abort ("Cannot solve generalized eigensystem") ;
    if (curv_rpt (repnum)) cout << "\n" ;

    spectrum = sqrt (spectrum) ;
    if (!has_nan (spectrum)) spectrum = sort (spectrum) ;
    lastrep  = repnum ;
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
    uint	blab	{ Blab::blablevel[BLAB::NUMERICS] } ;
    uint	offset	{ global.stage ? global.nobsG() : 0 } ;
    const auto&	geos	{ global.data().geos } ;
    const auto&	bckt	{ global.info().bckt } ;
    uint	nobs	{ global.info().nobs } ;
    int		ngens	{ global.info().neven.front() } ;

//    std::sort (bckt.begin(), bckt.end(), [geos](const uint3& a, const uint3& b)
//	{ return geos[a[0]].entry().filepos < geos[b[0]].entry().filepos ; }) ;

    if (global.interrupt) return ;
    if (ncol(numerics.delta) != ngens) gripe ("\nNeed to (re)evaluate delta!"); 
    if (blab > 2) cout << format("do_dvev: s = {:.6f}, nvev = {}, ",s,ncol(v)) ;

    set_zero (dv, nobs) ;
    numerics.dvev_buf = memptr(dv) ;

    if (memptr(v) != &numerics.vev[offset])
	{
	numerics.vev_tmp = numerics.vev ;
	if (global.stage)	numerics.vev_tmp.tail (ncol(v)) = v ;
	else			numerics.vev_tmp.head (ncol(v)) = v ;
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
	auto	nvev { ncol(v) } ;

	for (int indx(global.stage ? 0 : 1) ; indx < nvev ; ++indx)
	    {
	    if (abs(v[indx]) > maxv)
		{
		maxv = abs(v[indx]) ;
		maxi = indx ;
		}
	    if (abs(dv[indx]) > maxdv)
		{
		maxdv = abs(dv[indx]) ;
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
    uint	offset	{ global.stage ? global.nobsG() : 0 } ;
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

    uint	blab { Blab::blablevel[BLAB::NUMERICS] } ;
    auto&	eps  { numerics.odetol } ;
    ulong	maxi { index_max (abs (err)) } ;
    //doub	norm { l2norm  (err) / sqrt (ncol(err)) } ;
    //doub	norm { infnorm (err / (abs (err) + eps)) } ;
    doub	norm { infnorm (err) } ;
    if (blab > 2) cout << "err_norm: " << norm << "\n" << flush ;
    return norm ;
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

void Numerics::status_rpt (uint iters, uint steps)
    {
    auto&	svdcut	{ numerics.svdcut } ;
    auto&	okneg	{ global.oknegeig } ;
    auto&	maxi	{ status.maxloopi  } ;
    auto&	maxv	{ status.maxloopv  } ;
    doub&	mostneg { status.negeig[0] } ;
    bool	euclidF { global.stage && theory.euclid } ;
    bool	unphys  { !global.stage && check_loops() } ;

    if (!curv_rpt() && !unphys) cout << " ok" ;
    else if (unphys) cout << format(" unphys (#{}={:3.1f})", maxi, maxv) ;
    cout << " " << iters << "/" << steps << " iters/steps\n" << flush ;

    if (!euclidF && !okneg && mostneg < -svdcut) abort ("Negative curvature eigenvalues") ;
    }

void Numerics::write_data (doub value)				// Write data to MMAfile
    {
    ofstream& out { global.MMAstream } ;

    if (!global.MMAfile.size())
	 global.MMAfile = Global::dfltfilename("m") ;
    if (!out.is_open())
	{
	auto	mode	{ std::ios::out } ;
	string	okmsg	{ "Writing results to " } ;
	string	path	{ global.MMAfile } ;

	if (global.MMAappend)
	    {
	    mode |= std::ios::app ;
	    okmsg = "Appending results to " ;
	    }
	if (path.find ('/') == path.npos)
	    path = global.MMAdir + '/' + path ;

	out.open (path, mode) ;
	string msg { out.good() ? okmsg : "Cannot write to " } ;
	cout << msg << path << "\n" ;
	}
    if (!out.good()) gripe ("Cannot write MMA results file!") ;

    data_write (out, theory.euclid ? "F" : "H", value, eval_H()) ;

    uint beg { global.stage ? global.nobsG() : 1 } ;
    uint end { global.info().MMAlimit + (global.stage ? beg : 1) } ;
    if  (end > global.nobs())
	 end = global.nobs() ;

    for (uint i(beg) ; i < end ; ++i)
	{
	data_write (out, ObsList::obs(i).print(), value, vev[i]) ;
	}
    for (auto i : global.info().MMAlist)
	{
	if (i >= beg && i < end) continue ;
	data_write (out, ObsList::obs(i).print(), value, vev[i]) ;
	}

    if (!theory.euclid)
	{
	for (int i(0) ; i < Rep::list.size() ; ++i)
	    {
	    if (built_rep(i))
		{
		data_write (out, Rep::list[i].name, value, eval_spectra(i)) ;
		}
	    }
	}
    out << flush ;
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

bool Numerics::check_loops ()						// Loop vevs < 1?
    {
    doub	maxv (0.0) ;
    uint	maxi (0) ;
    auto	nvev { numerics.vev.size() } ;

    for (int indx(1) ; indx < nvev ; ++indx)
	{
	auto p { ObsList::obs.at(indx) } ;
	if (p->is_Loop() && abs(numerics.vev[indx]) > maxv)
	    {
	    maxv = abs(numerics.vev[indx]) ;
	    maxi = indx ;
	    }
	}
    status.maxloopv = numerics.vev[maxi] ;
    status.maxloopi = maxi ;
    return maxv > 1.0 ;
    }

bool Numerics::built_rep (int repnum)					// Is irrep built?
    {
    const auto&	curv	{ global.data().curv[repnum] } ;
    const auto&	lagr	{ global.data().lagr[repnum] } ;
    int		ngens	( global.info().gens[repnum].size() ) ;
    int		Hterms	( global.info().Hterms.size() ) ;

    return  lagr.entry().nrow == ngens && lagr.entry().ncol == ngens &&
	    curv.entry().nrow == ngens && curv.entry().ncol == ngens && 
	    curv.entry().nslice == Hterms && ngens > 0 ;
    }

void Numerics::data_write (ofstream& out, const string s, doub c, doub v)		// Save data value
    {
    out << s << "[" << c << "] = " << MMAform(v) << " ;\n" ;
    }

void Numerics::data_write (ofstream& out, const string s, doub c, const Rvec& vec)	// Save data vector
    {
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

void Numerics::data_write (ofstream& out, const string s, doub c, const Cvec& vec)	// Save data vector
    {
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
    string s { buf.str() } ;
    auto pos { s.find ("e") } ;
    if (pos != s.npos) { s.replace(pos, 1, "*^") ; }
    return s ;
    }

void Numerics::numericsinit ()					// Initialize expectation values
    {
    numerics.vev.resize (global.nobs()) ;
    if (global.stage == Global::Gauge)
	{
	set_zero (numerics.vev) ;
	numerics.vev[0] = 1.0 ;
	}
    else // global.stage == Global::Fermi
	{
	char8	mass { "mass" } ;
	auto	m    { Coupling::list[Coupling::indx (mass)].value } ;
	doub	condensate ;
	if (theory.euclid)	condensate = -1.0 / m ;
	else			condensate = m > 0 ? -0.5 : 0.5 ;
	set_zero (numerics.vev.tail (global.nobsF())) ;
	for (const auto& [indx_f,indx_g] : ObsList::obs.fermiinit)
	    {
	    numerics.vev[indx_f] = numerics.vev[indx_g] * condensate ;
	    }
	}
    if (!numerics.rk.nstage)
	 numerics.rk = RKdef::list.back() ;
    }
