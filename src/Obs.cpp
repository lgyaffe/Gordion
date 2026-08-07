#include "Obs.h"
#include "Global.h"
#include "Numerics.h"
#include "Op.h"
#include "Save.h"
#include "Symm.h"
#include "Gripe.h"
#include "Blab.h"
#include <regex>

Obs::Obs(const Op& a)						// Convert Op -> Obs
    :
    Str(a), corder(a.order), xorder(-1)
    {
    type = a.is_Fermion() ? ObsType::Fermion :
	   a.is_Eloop()   ? ObsType::Eloop : ObsType::Loop ;
    }

Obs::Obs(const string& s)					// Construct Obs
    :
    Str(s), type(obstype(s)), corder(-1), xorder(-1)
    {
    if (Blab::level(Blab::OBS) > 3)
	{
	cout << "Obs(string&) " << s << "\n" << flush ;
	}
    joinends() ;
    validate() ;
    }

Obs::Obs(const string& s, ObsType t, short cord, short xord)	// Construct Obs
    :
    Str(s), type(t), corder(cord), xorder(xord)
    {
    if (Blab::level(Blab::OBS) > 3)
	{
	cout << "Obs(string&,short,short) " << s << "\n" << flush ;
	}
    joinends() ;
    validate() ;
    }

void Obs::validate () const					// Test if valid Obs
    {
    if (Blab::level(Blab::OBS) > 1)
	{
	cout << "Obs::validate " << *this << "\n" << flush ;
	}
    if (size() == 0 && is_Loop()) return ;

    int	S	(0) ;
    int	fs	(0) ;
    int	Fs	(0) ;
    int	Es	(0) ;
    int	nf	(0) ;
    int	dim	(0) ;
    int derivs	(0) ;
    bool iseuc	{ theory.euclid } ;

    const char* err { nullptr } ;
    for (const auto c : *this)
	{
	if (islink(c) && dim <= axis(c))	dim = axis(c)+1 ;
	if (isferm(c) && nf  <= flav(c))	nf  = flav(c)/2+1 ;
	if (is_S(c))				++S ;
	if (isferm(c) && !isconj(c))		++fs ;
	if (isferm(c) &&  isconj(c))		++Fs ;
	if (isferm(c) && isderiv(c))		++derivs ;
	if (isE(c)  || isElink(c))		++Es ;
	if (isEE(c) || isEElink(c))		Es += 2 ;
	}
    if      (nf  > theory.nf)	 	err = "has excess fermion flavors" ;
    else if (dim > theory.dim)	 	err = "exceeds lattice dimension" ;
    else if (S && size() != 1)	 	err = "is malformed entropy" ;
    else if (fs + Fs && (fs != Fs))	err = "has bad fermion insertions" ;
    else if (fs + Fs > 2)		err = "has excessive fermions" ;
    else if (fs + Es > 2)		err = "has excessive E's" ;
    else if (fs == 0 && !isclosed())	err = "is not closed loop" ;
    else if (fs && !is_F(front()))	err = "is malformed fermion bilinear" ;
    else if (fs && !is_f(back()))	err = "is malformed fermion bilinear" ;
    else if (Es && islink(front()))	err = "is mis-rotated" ;
    else if (iseuc && derivs)		err = "has Grassmann derivative" ;
    else if (fs &&  Es  && type != ObsType::Efermion) err = "has wrong type" ;
    else if (fs && !Es  && type != ObsType::Fermion)  err = "has wrong type" ;
    else if (!fs && !Es && type != ObsType::Loop)     err = "has wrong type" ;
    else if (!fs && Es==1 && type != ObsType::Eloop)  err = "has wrong type" ;
    else if (!fs && Es>1  && type != ObsType::EEloop) err = "has wrong type" ;

    if (!err && size() > 1)
	{
	if (ligature (back(),front()))	err = "is malformed" ;
	for (auto p = cbegin()+1 ; !err && p < cend() ; ++p)
	    {
	    if (ligature (p[-1],p[0])) err = "is malformed" ;
	    }
	}
    if (err) gripe (format("Bad Obs: {} {}", print(), err)) ;
    }

short Obs::middleE() const				// Set & return midE location
    {
    if (midE < 0)
	{
	if (!doubleE(front()))
	    {
	    for (auto ptr = cbegin() ; ++ptr < cend() ;)
		{
		if (singleE (*ptr)) return midE = ptr - cbegin() ;
		}
	    fatal (format ("Missing midE in obs {}, type {}", print(), (int)type)) ;
	    }
	else midE = 0 ;
	}
    return midE ;
    }

bool Obs::Esublat() const				// E sublattice ?
    {
    short midE { middleE() } ;
    symb  c { (*this)[midE] } ;
    return is_Efermion() ? !(midE % 2) ^ (is_LE(c) || is_LEl(c)) : false ;
    }

ObsType Obs::obstype (const string s)			// Determine Obs type
    {
    std::regex			E_s ("[ABCD]") ;
    std::regex			F_s ("[FGHI]") ;
    std::sregex_iterator	end ;
    std::sregex_iterator	Ebeg (s.begin(), s.end(), E_s) ;
    std::sregex_iterator	Fbeg (s.begin(), s.end(), F_s) ;
    int				nE ( std::distance (Ebeg, end) ) ;
    int				nF ( std::distance (Fbeg, end) ) ;

    if      (nE == 0 && nF == 0) return ObsType::Loop ;
    else if (nE == 1 && nF == 0) return ObsType::Eloop ;
    else if (nE == 2 && nF == 0) return ObsType::EEloop ;
    else if (nE == 0 && nF == 1) return ObsType::Fermion ;
    else if (nE == 1 && nF == 1) return ObsType::Efermion ;
    else gripe (format ("Obs {} is invalid type", s)) ;
    }

int Obs::findstart() noexcept				// Rotate to preferred start
    {
    if (is_Loop() && size() > 0)
	{
	int a(0) ;
	int b(0) ;
	int len ( size() ) ;
	auto s { c_str() } ;

	while (++b < len)
	    {
	    int k(0) ;
	    while (k < len && s[(a+k) % len] == s[(b+k) % len]) ++k ;
	    if (k < len && s[(b+k) % len] < s[(a+k) % len]) a = b ;
	    }
	if (a) rotate (begin(), begin() + a, end()) ;
	}
    else if (is_EEloop() && !doubleE(front()))
	{
	const auto&	blab { Blab::level(Blab::OBS) } ;
	int		a(0), b(middleE()), k(0) ;
	int		len ( size() ) ;
	auto		s { c_str() } ;

	while (k < len && s[(a+k) % len] == s[(b+k) % len]) ++k ;
	if (k < len && s[(b+k) % len] < s[(a+k) % len])
	    {
	    if (blab > 2)
		{
		cout << "Obs::findstart: rotating obs " << *this << " "
		     << b << " places " << flush ;
		}
	    rotate (begin(), begin() + b, end()) ;
	    midE = len - b ;
	    if (blab > 2) cout << " -> " << *this << "\n" ;
	    }
	}
    return 1 ;
    }

int Obs::trans (const Symm& symm, int start) noexcept	// Symmetry transform Obs
    {
    if (front() == EntrG || front() == EntrF)	return 1 ;
    if (!start && symm.is_id())			return 1 ;

    const auto&	blab	{ Blab::level(Blab::SYMM) } ;
    bool	neg	{ false } ;
    if (blab > 1) cout << "Symm " << symm.name << " on " << *this
		       << " at " << start << " -> " ;
    if (symm.isCodd())
	{
	std::reverse (begin(), end()) ;
	start = size() - 1 - start ;
	midE  = -1 ;
	}
    if (start)
	{
	rotate (begin(), begin() + start, end()) ;
	midE = -1 ;
	}

    for (auto& c : *this)
	{
	neg ^= symm.sgn[c] ;
	c = symm.map[c] ;
	}
    if (!theory.euclid && isferm(back()) && isstag(back()))
	{
	neg ^= oddlen() ;
	back()  = stag(back()) ;
	front() = stag(front()) ;
	}
    if (blab > 1) cout << (neg ? "-" : "") << *this << "\n" << flush ;

    ++global.count().symmtrans ;
    return neg ? -1 : 1 ;
    }

ObsList::ObsList (const string s, bool can, bool clsfy)		// Construct ObsList
    : name(s), canonicalize(can), classify (clsfy)
    {
    store (Obs(Str(), ObsType::Loop, 0, 0)) ;
    }

PolyTerm ObsList::is_known (Obs&& a) const			// Find in ObsList
    {
    int		sgn  ( canonicalize ? a.canon() : a.findstart() ) ;
    numb	indx { find (a) } ;

    return (indx < MAXNUM-1) ? PolyTerm (PolyIndx(indx), sgn)
			     : PolyTerm (PolyIndx(),     0.0) ;
    }

PolyTerm ObsList::is_known (Obs&& a, Obs&& b) const		// Find in ObsList
    {
    int		sgna  ( canonicalize ? a.canon() : a.findstart() ) ;
    int		sgnb  ( canonicalize ? b.canon() : b.findstart() ) ;
    numb	indxa { find (a) } ;
    numb	indxb { find (b) } ;

    return (indxa < MAXNUM-1 && indxb < MAXNUM-1) ?
	PolyTerm (PolyIndx(indxa, indxb), sgna * sgnb) :
	PolyTerm (PolyIndx(), 0.0) ;
    }

PolyTerm ObsList::catalog (Obs a)			// Catalog Obs in list
    {							// N.B. pass by value
    const auto& blab { Blab::level(Blab::OBS) } ;
    if (blab > 1) cout << "catalog " << name << ": " << a << "\n" ;
    if (Obs::check) a.validate() ;

    int  sgn  { canonicalize ? a.canon() : a.findstart() } ;
    numb indx { find(a) } ;
    if (indx < MAXNUM-1)
	{
	if (blab > 1) cout << "catalog " << name << ": found "
			   << a << " at " << indx << "\n" ;
	return PolyTerm (indx, sgn) ;
	}
    else if (classify && !a.known_xord())
	{
	if (!a.classify (*this))
	    fatal (format("Failed classify: {} ({},{}) in {}",
		a.print(), a.corder, a.xorder, name)) ;
	}
    if (frozen())
	gripe (format("Cannot catalog {} into frozen list {}",
		a.print(), name)) ;
    PolyTerm ans (store(a), sgn) ;
    if (blab > 1) cout << "catalog " << name << ": stored " << a
		       << " -> " << ans[0] << "\n" ;
    return ans ;
    }

PolyTerm ObsList::catalog (Obs a, Obs b)		// Catalog Obs in list
    {							// N.B. pass by value
    const auto&	blab { Blab::level(Blab::OBS) } ;
    if (blab > 1) cout << "catalog " << name << ": " << a << ", " << b << "\n" ;
    if (Obs::check) { a.validate() ; b.validate() ; }

    int	sgna ( canonicalize ? a.canon() : a.findstart() ) ;
    int	sgnb ( canonicalize ? b.canon() : b.findstart() ) ;
    numb indxa { find(a) } ;
    numb indxb { find(b) } ;

    if (indxa < MAXNUM-1 && indxb < MAXNUM-1)
	{
	if (blab > 1)
	    cout << "catalog " << name << ": found "
		 << a << " at " << indxa << ", "
		 << b << " at " << indxb << "\n" ;
	return PolyTerm (PolyIndx (indxa,indxb), sgna * sgnb) ;
	}
    if (frozen())
	gripe (format("Cannot catalog {} {} into frozen list ()",
		    a.print(), b.print(), name)) ;
    if (classify && !a.known_xord())
	{
	if (!a.classify (*this))
	    fatal (format("Failed classify: {} ({},{}) in {}",
		a.print(), a.corder, a.xorder, name)) ;
	}
    if (classify && !b.known_xord())
	{
	if (!b.classify (*this))
	    fatal (format("Failed classify: {} ({},{}) in {}",
		b.print(), b.corder, b.xorder, name)) ;
	}
    PolyTerm ans (PolyIndx(store(a), store(b)), sgna * sgnb) ;
    if (blab > 1)
	{
	cout << "catalog " << name << ": stored/found " << a << " & "
	     << a << " -> " << ans[0] << " & " << ans[1] << "\n" ;
	}
    return ans ;
    }

numb ObsList::store (const Obs& o)			// Store Obs in ObsList
    {
    const auto&	blab { Blab::level(Blab::OBS) } ;
    if (classify)
	{
	if (o.bilinear() && !o.is_coord())
	    {
	    fatal (format("Non-coord Obs {} in {}", o.print(), name)) ;
	    }
	if (o.corder < 0)
	    {
	    cout << "Cannot store Obs with unknown corder\n" ;
	    }
	if (o.corder < 0 || o.xorder < 0)
	    {
	    fatal (format("Unclassified Obs {} ({},{}) in {}",
		o.print(), o.corder, o.xorder, name)) ;
	    }
	}
    if (blab > 1) cout << "Storing in " << name << ": " << o << "\n" ;
    auto [iter, isnew] { map.try_emplace (o, (*this).size()) } ;
    numb indx { iter->second } ;
    if (isnew)
	{
	push_back (&(iter->first)) ;
	if (!neq (ObsList::obs))
	    {
	    auto& info { global.info(o.is_fermi()) } ;
	    hasher (info.obshash, o) ;
	    if (info.nobs == MAXNUM-1)
		gripe ("Max # Obs exceeded: recompile without NUM32!") ;
	    ++info.nobs ;
	    }
	}
    return indx ;
    }

void ObsList::clear ()					// Empty list
    {
    if (!neq (ObsList::obs))
	{
	global.info(0).nobs = 0 ;
	global.info(1).nobs = 0 ;
	global.info(0).obshash = 0 ;
	global.info(1).obshash = 0 ;
	}
    Obsmap().swap (map) ;
    vector<const Obs*>().swap (*this) ;
    store (Obs(Str(), ObsType::Loop, 0, 0)) ;
    }

void ObsList::purge (numb limit)			// Purge entries
    {
    resize (limit) ;
    if (!neq (ObsList::obs))
	{
	auto& nobsG { global.info(0).nobs } ;
	auto& nobsF { global.info(1).nobs } ;
	nobsF -= std::erase_if (map, [&](const auto& p)
	    { return p.second >= limit && p.second >= nobsG; }) ;
	nobsG -= std::erase_if (map, [&](const auto& p)
	    { return p.second >= limit && p.second <  nobsG; }) ;
	if (nobsG + nobsF != size())
	    abort ("purge: Inconsistent ObsList size") ;
	rehash () ;
	}
    }

void ObsList::rehash () const			// Recalculate list hashes
    {
    global.info(0).obshash = 0 ;
    global.info(1).obshash = 0 ;
    for (const auto& ptr : *this)
	{
	hasher (global.info(ptr->is_fermi()).obshash, *ptr) ;
	}
    }

void ObsList::hasher (ulong& hash, const Obs& o) const	// List hasher
    {
    ulong x { std::hash<string>{}(o) } ;
    hash ^= x + 0x9e3779b9 + (hash << 6) + (hash >> 2) ;
    }

void ObsList::ondisk ()			// Leave ObsList::obs on disk
    {
    auto	stage	{ global.stage } ;
    const auto&	infoG	{ global.info(0) } ;
    const auto&	infoF	{ global.info(1) } ;
    auto&	list	{ ObsList::obs } ;
    auto	nobsG	{ infoG.nobs } ;
    auto	nobsF	{ infoF.nobs } ;

    if (list.size() > 1) cout << "Swapping master Obs list to disk\n" ;
    if (nobsG > 1)
	{
	bool		is_open	{ infoG.sysfile.stream.is_open() } ;
	const auto&	obs	{ global.data(0).obs } ;
	const auto&	nobs	{ obs.entry().items() } ;
	global.stage = Global::Gauge ;
	if (!is_open || nobsG != nobs) Save::save_sys () ;
	global.stage = stage ;
	}
    if (nobsF > 0)
	{
	bool		is_open	{ infoF.sysfile.stream.is_open() } ;
	const auto&	obs	{ global.data(1).obs } ;
	const auto&	nobs	{ obs.entry().items() } ;
	global.stage = Global::Fermi ;
	if (!is_open || nobsF != nobs) Save::save_sys () ;
	global.stage = stage ;
	}
    list.clear() ;
    ObsList::swapped = true ;
    }

void ObsList::obsinit (int stage)		// Load basic Obs
    {
    bool iseuc	 { theory.euclid } ;
    char link[4] { 'x', 'y', 'z', 'w' } ;
    char Link[4] { 'X', 'Y', 'Z', 'W' } ;

    freeze = false ;
    if (stage == 0)
	{
	if (iseuc)				// gauge entropy
	    {
	    Str entG ( EntrG ) ;
	    catalog (Obs(entG,ObsType::Entropy,0,0)) ;
	    }
	if (!iseuc)				// E & EE
	    {
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		Str E ( justE + i ) ;
		catalog (Obs(E,ObsType::Eloop,0,4)) ;
		Str EE ( KinG + i ) ;
		catalog (Obs(EE,ObsType::EEloop,0,4)) ;
		}
	    }
	for (int i(0) ; i < theory.dim ; ++i)	// single plaquettes
	    {
	    for (int j(i) ; ++j < theory.dim ;)
		{
		string plaq {link[i],link[j],Link[i],Link[j]} ;
		string Plaq {link[i],Link[j],Link[i],link[j]} ;
		catalog (Obs(plaq,ObsType::Loop,2,2)) ;
		catalog (Obs(Plaq,ObsType::Loop,2,2)) ;
		}
	    }
	for (int i(0) ; i < theory.dim ; ++i)	// Polyakov loops 
	    {
	    if (theory.box.comp[i])
		{
		string polyakov (theory.box.comp[i], char('x'+i)) ;
		string Polyakov (theory.box.comp[i], char('X'+i)) ;
		catalog (Obs(polyakov,ObsType::Loop,2,2)) ;
		catalog (Obs(Polyakov,ObsType::Loop,2,2)) ;
		}
	    }
	if (!neq (ObsList::obs)) global.info(0).maxord = 2 ;
	}
    else if (theory.nf)
	{
	if (iseuc)				// fermionic entropy
	    {
	    Str entF ( EntrF ) ;
	    catalog (Obs(entF,ObsType::Entropy,0,0)) ;
	    }
	if (!iseuc)				// gauge kinetic, fermion contrib
	    {
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		Str ee ( KinF + i ) ;
		catalog (Obs(ee,ObsType::EEloop,0,2)) ;
		}
	    }
	for (int i(0) ; i < theory.nf ; i+=2)	// fermion mass
	    {
	    string Ff {char('F'+i),char('f'+i)} ;
	    string Gf {char('G'+i),char('f'+i)} ;
	    if (iseuc)	catalog (Obs(Ff,ObsType::Fermion,0,0)) ;
	    else	catalog (Obs(Gf,ObsType::Fermion,0,0)) ;
	    }
	for (int i(0) ; i < theory.nf ; i+=2)	// fermion hopping
	    {
	    for (int j(0) ; j < theory.dim ; ++j)
		{
		string Fxf {char('F'+i),link[j],char('f'+i)} ;
		string FXf {char('F'+i),Link[j],char('f'+i)} ;
		catalog (Obs(Fxf,ObsType::Fermion,1,1)) ;
		catalog (Obs(FXf,ObsType::Fermion,1,1)) ;
		}
	    }
	if (!neq (ObsList::obs))
	    {
	    global.info(1).maxord = 1 ;
	    do_fermiinit() ;
	    }
	}
    freeze = true ;

    if (!neq (ObsList::obs) && global.info(0).nobs + global.info(1).nobs != size())
	abort (format("Inconsistent ObsList size: {} + {} != {}",
	    global.info(0).nobs, global.info(1).nobs, size())) ;
    }

int ObsList::do_fermiinit ()			// Initialize fermion -> loop map
    {
    const auto&	blab	{ Blab::level(Blab::OBS) } ;
    long	beg	{ global.info(0).nobs } ;
    uint	initfail (0) ;
    if (blab > 3) cout << "do_fermiinit start\n" << flush ;

    fermiinit.clear () ;
    for (long i(beg) ; i < size() ; ++i)
	{
	const Obs& a { (*this)(i) } ;
	if (!a.is_Fermion() || !a.isclosed()) continue ;

	Obs b (a.cbegin()+1, a.cend()-1, ObsType::Loop) ;
	b.joinends() ;
	if (blab > 3) cout << "fermion " << a << " loop partner " << b << "\n" ;
	PolyTerm term { is_known (std::move(b)) } ;
	if (!term.coeff)
	    {
	    ++initfail ;
	    if (blab > 1) cout << "Warning: Cannot initialize " << a << "\n" ;
	    }
	else fermiinit.emplace_back (i, term[0]) ;
	}
    return initfail ;
    }

ostream& ObsList::print (ostream& stream, numb indx) const	// Print indexed Obs
    {
    bool	addvev	 { !neq(ObsList::obs) } ;
    auto	prevprec { stream.precision(12) } ;
    const Obs&	obs	 { *at(indx) } ;
    stream << name << " Obs #" << indx << ": " << std::setprecision(12) ;
    stream << obs ;
    if (addvev && indx < numerics.vev.size())
	stream << " = " << numerics.vev[indx] ;
    stream << "\n" << std::setprecision (prevprec) ;
    return stream ;
    }

ostream& operator<< (ostream& stream, const Obs& obs)		// Print Obs -> stream
    {
    if (Blab::level(Blab::OBS))
	{
	stream  << "(" << obs.corder << "," << obs.xorder ;
	stream << "," << Obs::type_name [(int)obs.type] ;
	stream << ") " ;
	}
    if (obs.imag()) stream << "i " ;
    stream << (Str) obs ;
    return stream ;
    }

ostream& ObsList::print (ostream& stream) const			// Print ObsList -> stream
    {
    bool	addvev	 { !neq(ObsList::obs) } ;
    auto	prevprec { stream.precision(12) } ;
    stream << name << " observables:\n" << std::setprecision(12) ;
    for (int indx(0) ; indx < size() ; ++indx)
	{
	stream << " #" << indx << ": " << *at(indx) ;
	if (addvev && indx < numerics.vev.size())
	    stream << " = " << numerics.vev[indx] ;
	stream << "\n" ;
	}
    stream << std::setprecision (prevprec) ;
    return stream ;
    }

ObsStats::ObsStats (const ObsList& list)			 // Construct ObsStats
    {
    maxloop = 0 ;
    doub	maxloopvev { 0.0 } ;
    ulong	lengthsum (0) ;
    auto&	vevs { numerics.vev } ;
    auto	nvev { vevs.size() } ;
    auto 	beg  { list.begin() } ;
    auto 	end  { list.end()   } ;

    for (auto ptr { beg } ; ptr < end ; ++ptr)
	{
	int	indx	( ptr - beg ) ;
	auto	obs	{ list(indx) } ;
	int	t	{ (int) obs.type } ;
	short	c	( std::max (obs.corder, (short) 0) ) ;
	short	x	( std::max (obs.xorder, (short) 0) ) ;
	
	if (c >= (*this)[t].size())	(*this)[t].resize (c+1) ;
	if (x >= (*this)[t][c].size())	(*this)[t][c].resize (x+1) ;

	if (c > (*this).maxc) (*this).maxc = c ;
	if (x > (*this).maxx) (*this).maxx = x ;

	++(*this)[t][c][x] ;
	lengthsum += obs.size() ;
	if (obs.is_Loop() && indx && indx < nvev)
	    {
	    if (std::abs(vevs[indx]) > maxloopvev)
		{
		maxloopvev = std::abs(vevs[indx]) ;
		maxloop = indx ;
		}
	    }
	}
    (*this).avglen = (float) lengthsum / list.size() ;
    }
