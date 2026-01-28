#include "Obs.h"
#include "Op.h"
#include "Symm.h"
#include "Global.h"
#include "Numerics.h"
#include "Gripe.h"
#include "Blab.h"
#include <regex>

Obs::Obs(const Op& a)						// Convert Op -> Obs
    : SymbStr(a), corder(a.order), xorder(-1)
    {
    type = a.is_Fermion() ? ObsType::Fermion :
	   a.is_Eloop()   ? ObsType::Eloop : ObsType::Loop ;
    }

Obs::Obs(const string& s)					// Construct Obs
    : SymbStr(to_symb(s)), type(obstype(s)), corder(-1), xorder(-1)
    {
    if (Blab::blablevel[BLAB::OBS] > 3)
	{
	cout << "Obs(string&) " << s << "\n" << flush ;
	}
    validate() ;
    }

Obs::Obs(const string&  s, ObsType t, short cord, short xord)	// Construct Obs
    : SymbStr(to_symb(s)), type(t), corder(cord), xorder(xord)
    {
    if (Blab::blablevel[BLAB::OBS] > 3)
	{
	cout << "Obs(string&,short,short) " << s << "\n" << flush ;
	}
    validate() ;
    }

void Obs::validate () const					// Test if valid Obs
    {
    if (Blab::blablevel[BLAB::OBS] > 1)
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
	if (isentpy(c))				++S ;
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
    else if (fs && !isF(front()))	err = "is malformed fermion bilinear" ;
    else if (fs && !isf(back()))	err = "is malformed fermion bilinear" ;
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
	if (!EEorEElink(front()))
	    {
	    for (auto ptr = cbegin() ; ++ptr < cend() ;)
		{
		if (EorElink(*ptr)) return midE = ptr - cbegin() ;
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
    return is_Efermion() ? !(midE % 2) ^ (isLE(c) || isLEl(c)) : false ;
    }

ObsType Obs::obstype (const string s)			// Determine Obs type
    {
    int				nE = 0 ;
    int				nF = 0 ;
    std::regex			E_s ("[ABCD]") ;
    std::regex			F_s ("[FGHI]") ;
    std::sregex_iterator	Ebeg (s.begin(), s.end(), E_s) ;
    std::sregex_iterator	Fbeg (s.begin(), s.end(), F_s) ;
    std::sregex_iterator	Eend ;
    std::sregex_iterator	Fend ;

    while (Ebeg != Eend) { ++nE ; ++Ebeg ; }	// count E's
    while (Fbeg != Fend) { ++nF ; ++Fbeg ; }	// count F's

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
    else if (is_EEloop() && !EEorEElink(front()))
	{
	uint	blab { Blab::blablevel[BLAB::OBS] } ;
	int	a(0), b(middleE()), k(0) ;
	int	len ( size() ) ;
	auto	s { c_str() } ;

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

    bool neg	{ false } ;
    uint blab	{ Blab::blablevel[BLAB::SYMM] } ;
    if (blab > 0) cout << "Symm " << symm.name << " on " << *this
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
    if (blab > 0) cout << (neg ? "-" : "") << *this << "\n" << flush ;

    ++global.count().symmtrans ;
    return neg ? -1 : 1 ;
    }

ObsList::ObsList (const string s, bool can, bool clas)		// Construct ObsList
    : name(s), canonicalize(can), classify (clas)
    {
    store (Obs(SymbStr(), ObsType::Loop, 0, 0)) ;
    }

PolyTerm ObsList::is_known (Obs&& a) const			// Find in ObsList
    {
    doub	sgn  ( canonicalize ? a.canon() : a.findstart() ) ;
    uint	indx { find (a) } ;

    return (indx < UINT_MAX-1) ? PolyTerm (Polyindx(indx), sgn)
			       : PolyTerm (Polyindx(),     0.0) ;
    }

PolyTerm ObsList::is_known (Obs&& a, Obs&& b) const		// Find in ObsList
    {
    doub	sgna  ( canonicalize ? a.canon() : a.findstart() ) ;
    doub	sgnb  ( canonicalize ? b.canon() : b.findstart() ) ;
    uint	indxa { find (a) } ;
    uint	indxb { find (b) } ;

    return (indxa < UINT_MAX-1 && indxb < UINT_MAX-1) ?
	PolyTerm (Polyindx(indxa, indxb), sgna * sgnb) :
	PolyTerm (Polyindx(), 0.0) ;
    }

PolyTerm ObsList::catalog (Obs a)			// Catalog Obs in list
    {							// N.B. pass by value
    uint blab { Blab::blablevel[BLAB::OBS] } ;
    if (blab > 1) cout << "catalog " << name << ": " << a << "\n" ;
    if (Obs::check) a.validate() ;

    int  sgn  { canonicalize ? a.canon() : a.findstart() } ;
    uint indx { find(a) } ;
    if (find (a) < UINT_MAX-1)
	{
	if (blab > 1)
	    cout << "catalog " << name << ": found "
		 << a << " at " << indx << "\n" ;
	return PolyTerm (indx, sgn) ;
	}
    else if (classify && !a.known_xord())
	{
	if (!a.classify (*this))
	    fatal (format("Failed classify: {} ({},{}) in {}",
		a.print(), a.corder, a.xorder, name)) ;
	}
    PolyTerm ans (store(a), sgn) ;
    if (blab > 1)
	{
	cout << "catalog " << name << ": stored " << a
	     << " -> " << ans[0] << "\n" ;
	}
    return ans ;
    }

PolyTerm ObsList::catalog (Obs a, Obs b)		// Catalog Obs in list
    {							// N.B. pass by value
    uint blab { Blab::blablevel[BLAB::OBS] } ;
    if (blab > 1) cout << "catalog " << name << ": " << a << ", " << b << "\n" ;

    int	sgna ( canonicalize ? a.canon() : a.findstart() ) ;
    int	sgnb ( canonicalize ? b.canon() : b.findstart() ) ;

    if (Obs::check) { a.validate() ; b.validate() ; }
    if (classify)
	{
	if (!a.known_xord()) a.classify (*this) ;
	if (!b.known_xord()) b.classify (*this) ;
	}
    PolyTerm ans (Polyindx(store(a), store(b)), sgna * sgnb) ;
    if (blab > 1)
	{
	cout << "catalog " << name << ": stored/found " << a << " & "
	     << a << " -> " << ans[0] << " & " << ans[1] << "\n" ;
	}
    return ans ;
    }

uint ObsList::store (const Obs& o)			// Store Obs in ObsList
    {
    if (classify)
	{
	if (o.corder < 0)
	    {
	    cout << "Cannot store Obs with unknown corder\n" ;
	    }
	if (o.bilinear() && !o.is_coord())
	    {
	    fatal (format("Non-coord Obs {} in {}", o.print(), name)) ;
	    }
	if (o.corder < 0 || o.xorder < 0)
	    {
	    fatal (format("Unclassified Obs {} ({},{}) in {}",
		o.print(), o.corder, o.xorder, name)) ;
	    }
	}
    auto [iter, isnew] { map.try_emplace (o, (*this).size()) } ;
    uint indx { iter->second } ;
    if (isnew) push_back (&(iter->first)) ;
    return indx ;
    }

void ObsList::obsinit ()			// Load basic Obs
    {
    bool iseuc	 { theory.euclid } ;
    int  stage	 { global.stage } ;
    char link[4] { 'x', 'y', 'z', 'w' } ;
    char Link[4] { 'X', 'Y', 'Z', 'W' } ;

    ObsList::freeze = false ;
    if (stage == 0)
	{
	if (iseuc)				// gauge entropy
	    {
	    SymbStr entG ( EntrG ) ;
	    catalog (Obs(entG,ObsType::Entropy,0,0)) ;
	    }
	if (!iseuc)			// E & EE
	    {
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		SymbStr E ( justE + i ) ;
		catalog (Obs(E,ObsType::Eloop,0,4)) ;
		SymbStr EE ( KinG + i ) ;
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
	if (!neq(ObsList::obs)) global.nobsG() = nobs() ;
	}
    else if (theory.nf)
	{
	if (iseuc)				// fermionic entropy
	    {
	    SymbStr entF ( EntrF ) ;
	    catalog (Obs(entF,ObsType::Entropy,0,0)) ;
	    }
	if (!iseuc)				// gauge kinetic, fermion contrib
	    {
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		SymbStr ee ( KinF + i ) ;
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
	if (!neq(ObsList::obs))
	    {
	    global.nobsF() = ObsList::obs.size() - global.nobsG() ;
	    do_fermi_init() ;
	    }
	}
    ObsList::freeze = true ;
    }

void ObsList::do_fermi_init ()			// Initialize fermion -> loop map
    {
    uint	initfail  ( 0 ) ;
    uint	beg	  { global.nobsG() } ;
    uint	blab	  { Blab::blablevel[BLAB::OBS] } ;
    if (blab > 3) cout << "do_fermi_init start\n" << flush ;

    for (uint i(beg) ; i < nobs() ; ++i)
	{
	const Obs& a { (*this)(i) } ;
	if (!a.is_Fermion() || !a.isclosed()) continue ;

	Obs b (a.cbegin()+1, a.cend()-1, ObsType::Loop) ;
	b.joinends() ;
	if (blab > 3) cout << "fermion " << a << " loop partner " << b << "\n" ;
	auto term { is_known (std::move(b)) } ;
	if (!term.coeff)
	    {
	    ++initfail ;
	    if (blab > 1) cout << "Warning: Cannot initialize " << a << "\n" ;
	    }
	else fermiinit.emplace_back (i, term[0]) ;
	}
    if (blab)
	{
	cout << "Fermion initializations: " << fermiinit.size() ;
	if (initfail) cout << " + " << initfail << " missing loop partners" ;
	cout << "\n" ;
	}
    }

ostream& ObsList::print (ostream& stream, uint indx) const	// Print indexed Obs
    {
    bool	addvev	{ !neq(ObsList::obs) } ;
    const Obs&	obs	{ *at(indx) } ;
    stream << name << " Obs #" << indx << ": " ;
    stream << obs ;
    if (addvev && indx < numerics.vev.size())
	stream << " = " << numerics.vev[indx] ;
    stream << "\n" ;
    return stream ;
    }

ostream& operator<< (ostream& stream, const Obs& obs)		// Print Obs -> stream
    {
    if (Blab::blablevel[BLAB::OBS])
	{
	stream  << "(" << obs.corder << "," << obs.xorder ;
	//if (obs.midE >= 0) stream << "," << obs.midE ;
	stream << "," << Obs::type_name [(int)obs.type] ;
	stream << ") " ;
	}
    if (obs.imag()) stream << "i " ;
    stream << (SymbStr) obs ;
    return stream ;
    }

ostream& ObsList::print (ostream& stream) const			// Print ObsList -> stream
    {
    bool addvev { !neq(ObsList::obs) } ;
    cout << name << " observables:\n" ;
    for (int indx(0) ; indx < size() ; ++indx)
	{
	stream << " #" << indx << ": " << *at(indx) ;
	if (addvev && indx < numerics.vev.size())
	    stream << " = " << numerics.vev[indx] ;
	stream << "\n" ;
	}
    return stream ;
    }

ObsStats::ObsStats (const ObsList& list)			 // Construct ObsStats
    {
    maxloop = 0 ;
    doub	maxloopvev { 0.0 } ;
    ulong	lengthsum (0) ;
    auto&	vevs { numerics.vev } ;
    auto	nvev { vevs.size() } ;
    auto	nobs { list.size() } ;
    for (int indx(0) ; indx < nobs ; ++indx)
	{
	auto	p { list.at(indx) } ;
	int	t { (int) p->type } ;
	short	c ( std::max (p->corder, (short) 0) ) ;
	short	x ( std::max (p->xorder, (short) 0) ) ;
	
	if (c >= (*this)[t].size())	(*this)[t].resize (c+1) ;
	if (x >= (*this)[t][c].size())	(*this)[t][c].resize (x+1) ;

	if (c > (*this).maxc) (*this).maxc = c ;
	if (x > (*this).maxx) (*this).maxx = x ;

	++(*this)[t][c][x] ;
	lengthsum += p->size() ;
	if (p->is_Loop() && indx && indx < nvev)
	    {
	    if (abs(vevs[indx]) > maxloopvev)
		{
		maxloopvev = abs(vevs[indx]) ;
		maxloop = indx ;
		}
	    }
	}
    (*this).avglen = (float) lengthsum / list.size() ;
    }
