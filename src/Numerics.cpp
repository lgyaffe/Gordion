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
    int  delta ;
    doub tol   (0) ;
    uint iters (0) ;
    uint steps (0) ;
    ++stats.tries ;
    status.reset() ;

    do  {
	if (global.interrupt || ++iters > minmax) break ;
	delta = do_step (tol) ;
	steps += abs(delta) ;
	tol = mintol ;
	} while (delta > 0) ;

    status_rpt (iters,steps) ;
    if (iters > minmax)
	{
	++stats.fails ;
	abort ("Minimization failed, max iterations exceeded") ;
	}
    else if (delta < 0)
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
    doub	s (0) ;
    bool	ok ;

    eval_delta() ;
    if (blab > 1) cout << "  |delta| = " << arma::norm (delta,2) << ", " ;
    if (blab > 2) cout << "\n   delta = \n" << delta << flush ;

    if (tol && arma::norm (delta,2) <= tol) return 0 ; // converged?

    if (global.info().nobs != global.nobs())		// integrate vev subvec
	{
	real*	vevptr { &vev [global.stage ? global.nobsG() : 0] } ; 
	Numvec	subvev { vevptr, nvev, false, true } ;
	ok = ode.integrate (s, 1.0, subvev) ;
	}
    else ok = ode.integrate (s, 1.0, vev) ;

    if (blab > 1) cout << ode.steps << " step(s) " << ode.rejects << " rejects "
		       << (theory.euclid ? "F = " : "H = ") << eval_H() << "\n" << flush ;

    return ok ? ode.steps : -ode.steps ;
    }

Numvec& Numerics::eval_delta (bool print)		// Evaluate delta vector
    {
    auto	opts { arma::solve_opts::refine + arma::solve_opts::equilibrate } ;
    Numvec&	grad { eval_grad() } ;
    Nummtx&	curv { eval_curv(0, false) } ;

    if (minlim && minlim < global.info().maxgen)
	{
	const auto&	gens { global.info().gens.front() } ;
	uvec		use  ( curv.n_rows ) ;
	int		n(0) ;

	for (int i(0) ; i < curv.n_rows && i < gens.size() ; ++i)
	    {
	    if (gens[i].order <= minlim) use(n++) = i ;
	    }
	use.resize (n) ;
	check_curv (curv(use,use)) ;
	Numvec del ;

	if (svdlim)	del = -arma::pinv  (curv(use,use), svdlim) * grad(use) ;
	else		del = -arma::solve (curv(use,use), grad(use), opts) ;
	delta.zeros (curv.n_rows) ;
	delta(use) = del ;
	}
    else
	{
	check_curv (curv) ;
	if (svdlim)	delta = -arma::pinv  (curv, svdlim) * grad ;
	else		delta = -arma::solve (curv, grad, opts) ;
	}
    if (print) cout << "Delta = \n" << delta ;
    return delta ;
    }

Numvec& Numerics::eval_grad (bool print)		// Evaluate gradient vector
    {
    const auto& coup	{ Coupling::list } ;
    const auto& grad	{ global.data().grad } ;
    const auto& Hterms	{ global.info().Hterms } ;
    int		ngens	{ global.info().neven.front() } ;
    int		nterms	( Hterms.size() ) ;
    int		Hterm   { -1 } ;
    doub	coeff ;

    if (grad.entry().ncol != nterms ||
	grad.entry().nrow != ngens) gripe ("Need to (re)build gradient!") ;

    gradient.zeros (ngens) ;
    for (const auto& poly : grad)
	{
	const GradHdr&	info ( poly ) ;
	doub		val  ( 0.0 ) ;
	for (const auto& term : poly) val += termvalue (term) ;
	if (info.Hterm != Hterm)
	    {
	    Hterm = info.Hterm ;
	    int	exp  { Hterms[Hterm].exponent } ;
	    int	indx { Hterms[Hterm].coupindx } ;
	    coeff = exp ? std::pow (coup[indx].value, exp) : 1.0 ;
	    }
	gradient [info.gen] += coeff * val ;
	}
    if (print) cout << "Gradient = \n" << gradient ;
    return gradient ;
    }

Nummtx& Numerics::eval_curv (int repnum, bool all, int print)	// Evaluate curvature
    {
    const auto& coup	{ Coupling::list } ;
    const auto& curv	{ global.data().curv[repnum] } ;
    const auto& Hterms	{ global.info().Hterms } ;
    int		ngen	( global.info().gens[repnum].size() ) ;
    int		neven	{ global.info().neven[repnum] } ;
    int		nuse	( all ? ngen : neven ) ;
    int		nterms	( Hterms.size() ) ;
    int		Hterm   { -1 } ;
    doub	coeff ;

    if (curv.entry().nslice != nterms ||
	curv.entry().ncol != ngen  ||
	curv.entry().nrow != ngen) gripe ("Need to (re)build curvature!") ;

    curvature.zeros (nuse,nuse) ;
    for (const auto& poly : curv)
	{
	const CurvHdr&	info ( poly ) ;
	if (info.gen1 < nuse && info.gen2 < nuse)
	    {
	    doub	val  ( 0.0 ) ;
	    for (const auto& term : poly) val += termvalue (term) ;
	    if (info.Hterm != Hterm)
		{
		Hterm = info.Hterm ;
		int	exp  { Hterms[Hterm].exponent } ;
		int	indx { Hterms[Hterm].coupindx } ;
		coeff = exp ? std::pow (coup[indx].value, exp) : 1.0 ;
		}
	    curvature (info.gen1,info.gen2) += coeff * val ;
	    }
	}
    if (print > 1)
	{
	cout << Rep::list[repnum].name << " Curvature = \n" << curvature ;
	}
    if (print)
	{
	cout << Rep::list[repnum].name << " Curvature eigenvalues = \n"
	     << arma::sort (arma::eig_gen(curvature.submat(0,0,neven-1,neven-1))) ;
	if (ngen > neven)
	     cout << arma::sort (arma::eig_gen(curvature.submat(neven,neven,ngen-1,ngen-1))) ;
	}
    return curvature ;
    }

Nummtx& Numerics::eval_lagr (int repnum, bool print)	// Evaluate Lagrange bracket matrix
    {
    const auto& lagr	{ global.data().lagr[repnum] } ;
    int		ngens	( global.info().gens[repnum].size() ) ;

    if (lagr.entry().ncol != ngens ||
	lagr.entry().nrow != ngens) gripe ("Need to build Lagrange matrix!") ;

    lagrange.zeros (ngens,ngens) ;
    for (const auto& poly : lagr)
	{
	const LagrHdr&	info ( poly ) ;
	doub		val  ( 0.0 ) ;
	for (const auto& term : poly) val += termvalue (term) ;
	lagrange (info.gen1,info.gen2) = val ;
	}
    if (print) cout << Rep::list[repnum].name << " Lagrange bracket = \n" << lagrange ;
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
	doub	polyval (0.0) ;
	int	indx  { Hterm.coupindx } ;
	int	exp   { Hterm.exponent } ;
	doub	coeff { exp ? std::pow (coup[indx].value, exp) : 1.0 } ;

	for (const auto& polyterm : Hterm.cpoly)
	    {
	    polyval += termvalue (polyterm) ;
	    }
	H += coeff * polyval ;
	}
    if (print)
	{
	auto	prev { cout.precision(12) } ;
	cout << label << " = " << H << std::setprecision (prev) << "\n" ;
	}
    return H ;
    }

Cmplxv& Numerics::eval_spectra (int repnum, bool print)	// Evaluate particle spectrum
    {
    std::complex i	{ 0.0, 1.0 } ;
    Nummtx	 curv	{ eval_curv (repnum, true) } ;
    Nummtx&	 lagr	{ eval_lagr (repnum) } ;
    uvec	 use	( curv.n_rows ) ;

    if (global.symcurv)
	{
	curv += arma::trans(curv) ;
	curv /= 2.0 ;
	}
    constexpr doub tiny { 1.e-100 } ;
    constexpr doub gunk { 1.e-10 } ;
    constexpr doub huge { 1.e80 } ;
    lagr.diag() += tiny ;
    status.reset() ;

    if (speclim && speclim < global.info().maxgen)
	{
	const auto&	gens { global.info().gens[repnum] } ;
	auto		ptr  { use.begin() } ;

	for (int i(0) ; i < curv.n_rows && i < gens.size() ; ++i)
	    {
	    if (gens[i].order <= speclim) *ptr++ = i ;
	    }
	use.resize (ptr - use.begin()) ;
	check_curv (curv(use,use)) ;
	if (!arma::eig_pair (spectrum,modes,curv(use,use),lagr(use,use)))
	    {
	    abort ("Cannot solve generalized eigensystem") ;
	    }
	}
    else
	{
	check_curv (curv) ;
	if (!arma::eig_pair (spectrum,modes,curv,lagr))
	    {
	    abort ("Cannot solve generalized eigensystem") ;
	    }
	}
    if (curv_rpt (repnum)) cout << "\n" ;

    spectrum *= -i ;
    if (spectrum.n_elem) spectrum.clean (10 * arma::datum::eps * abs(spectrum.max())) ;

    use.resize (spectrum.n_elem) ;
    auto ptr { use.begin() } ;
    for (int i(0) ; i < spectrum.n_elem ; ++i)
	{
	auto eig { spectrum[i] } ;
	if (eig.real() < 0 || (eig.real() == 0 && eig.imag() < -gunk)) continue ;
	if (std::isfinite(abs(eig)) && abs(eig) < huge) *ptr++ = i ;
	}
    use.resize (ptr - use.begin()) ;
    spectrum = spectrum(use) ;
    modes    = modes.cols(use) ;
    lastrep  = repnum ;

    std::sort (spectrum.begin(), spectrum.end(), [](const cmplx& a, const cmplx& b)
	{
	if (abs(a.imag()) < gunk && abs(b.imag()) < gunk) return a.real() < b.real() ;
	else if (a.imag() < gunk)			  return true ;
	else if (b.imag() < gunk)			  return false ;
	else return abs(a.imag()) < abs(b.imag()) ;
	}) ;

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

void Numerics::do_dvev (doub s, const Numvec& v, Numvec& dv)	// Evaluate vev derivs
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
    if (numerics.delta.n_elem != ngens) gripe ("\nNeed to (re)evaluate delta!"); 
    if (blab > 2) cout << format("do_dvev: s = {:.6f}, ", s) ;

    dv.zeros (nobs) ;
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
	uint	maxi  (0) ;
	uint	maxdi (0) ;
	auto	nvev { numerics.vev.size() } ;

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
	cout << format("maxv #{}={:.2g}, ", maxi,  maxv)  ;
	cout << format("maxdv #{}={:.2g}\n", maxdi, maxdv) << flush ;
	}
    }

void Numerics::do_dvev_bckt (const uint3& bucket)		// Evaluate dvev bucket
    {
    if (global.interrupt) return ;

    uint	bcktnum	{ bucket[0] } ;
    uint	first	{ bucket[1] } ;
    uint	last	{ bucket[2] } ;
    uint	offset	{ global.stage ? global.nobsG() : 0 } ;
    auto&	geos	{ global.data().geos[bcktnum] } ;
    auto&	delta	{ numerics.delta } ;
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

	if (info.indx < first || info.indx > last)
	    fatal (format ("Inconsistent geo record: bucket {} indx {} first {} last {}!",
		    bcktnum, info.indx, first, last)) ; 

	if (global.interrupt) break ;
	for (const auto& term : poly) val += termvalue (term,v) ;
	dv [info.indx - offset] += delta[info.gen] * val ;
	}
    if (global.geoswap) Save::read_geo_bckt (-bcktnum-1) ;
    }

doub Numerics::err_norm (const Numvec& err, const Numvec& y)	// ODE error vector norm
    //
    //	Evaluate ODE error norm
    //
    {
    if (global.interrupt) return 0 ;

    uint	blab { Blab::blablevel[BLAB::NUMERICS] } ;
    doub	norm { arma::norm (err,2) / sqrt (err.n_elem) } ;
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
    auto&	svdlim	{ numerics.svdlim } ;
    auto&	okneg	{ global.oknegeig } ;
    auto&	maxi	{ status.maxloopi  } ;
    auto&	maxv	{ status.maxloopv  } ;
    doub&	mostneg { status.negeig[0] } ;
    bool	euclidF { global.stage && theory.euclid } ;
    bool	unphys  { !global.stage && check_loops() } ;

    if (!curv_rpt() && !unphys) cout << " ok" ;
    else if (unphys) cout << format(" unphys (#{}={:3.1f})", maxi, maxv) ;
    cout << " " << iters << "/" << steps << " iters/steps\n" << flush ;

    if (!euclidF && !okneg && mostneg < -svdlim) abort ("Negative curvature eigenvalues") ;
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

bool Numerics::check_curv (const Nummtx& curv)				// OK curvature eigs?
    {
    if (curv.n_rows) 
	{
	Cmplxv	eigs  { arma::eig_gen(curv) } ;
	uvec	reals { arma::sort_index (arma::real (eigs)) } ;
	uvec	imags { arma::sort_index (arma::imag (eigs),"descend") } ;
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

void Numerics::data_write (ofstream& out, const string s, doub c, const Numvec& vec)	// Save data vector
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

void Numerics::data_write (ofstream& out, const string s, doub c, const Cmplxv& vec)	// Save data vector
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
	numerics.vev.zeros() ;
	numerics.vev[0] = 1.0 ;
	}
    else // global.stage == Global::Fermi
	{
	doub	condensate (-0.5) ;
	if (theory.euclid)
	    {
	    char8 mass  { "mass" } ;
	    int   mindx	{ Coupling::indx (mass) } ;
	    condensate = -1.0 / Coupling::list[mindx].value ;
	    }
	numerics.vev.tail(global.nobsF()).zeros() ;
	for (const auto& [indx_f,indx_g] : ObsList::obs.fermiinit)
	    {
	    numerics.vev[indx_f] = numerics.vev[indx_g] * condensate ;
	    }
	}
    if (!numerics.rk.nstage)
	 numerics.rk = RKdef::list.back() ;
    }
