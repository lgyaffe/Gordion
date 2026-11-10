#include "List.h"
#include "Commute.h"
#include "Global.h"
#include "Numerics.h"
#include "Poly.h"
#include "Save.h"
#include "Blab.h"

ostream& OpList::print (ostream& stream, numb indx) const	// Print indexed Op
    {
    const Op& op { (*this)[indx] } ;
    return stream << " op #" << indx << " = " << op << "\n" ;
    }

ostream& OpList::print (ostream& stream) const			// Print OpList
    {
    stream << " operators:\n" ;
    for (int indx(0) ; indx < size() ; ++indx)
	{
	const Op& op { (*this)[indx] } ;
	stream  << " #" << indx << " = " << op << "\n" ;
	}
    return stream ;
    }

void OpList::opinit (uint stage)		// Operator initialization
    {
    bool isham { !theory.euclid } ;
    char l[4]  { 'x', 'y', 'z', 'w' } ;
    char L[4]  { 'X', 'Y', 'Z', 'W' } ;
    char f[4]  { 'f', 'g', 'h', 'i' } ;
    char F[4]  { 'F', 'G', 'H', 'I' } ;

    if (stage == 0)
	{
	OpType	loop { OpType::Loop } ;
	OpList&	list { global.info(0).ops } ;

	for (int i(0) ; i < theory.dim ; ++i)		// 1x1 plaq
	    {
	    for (int j(0) ; j < theory.dim ; ++j)
		{
		if (i == j) continue ;
		Op plaq { string {l[i],l[j],L[i],L[j]}, loop, 2 } ;
		if (isham) list.store (plaq) ;
		OpSum::loop_dt (plaq, list) ;

		Op Plaq { string {l[i],L[j],L[i],l[j]}, loop, 2 } ;
		if (isham) list.store (Plaq) ;
		OpSum::loop_dt (Plaq, list) ;
		}
	    }
	for (int i(0) ; i < theory.dim ; ++i)		// 2x1 fig-8, possibly bent
	    {
	    for (int j(0) ;  j < theory.dim ; ++j)
		{
		for (int k(0) ;  k < theory.dim ; ++k)
		    {
		    if (i == k || j == k) continue ;
		    Op figa { string {l[i],l[k],l[j],L[k],L[j],l[k],L[i],L[k]}, loop, 4 } ;
		    if (isham) list.store (figa) ;
		    OpSum::loop_dt (figa, list) ; 
		    Op figb { string {l[i],L[k],l[j],l[k],L[j],L[k],L[i],l[k]}, loop, 4 } ;
		    if (isham) list.store (figb) ;
		    OpSum::loop_dt (figb, list) ; 
		    Op figc { string {l[i],l[k],L[j],L[k],l[j],l[k],L[i],L[k]}, loop, 4 } ;
		    if (isham) list.store (figc) ;
		    OpSum::loop_dt (figc, list) ; 
		    Op figd { string {l[i],L[k],L[j],l[k],l[j],L[k],L[i],l[k]}, loop, 4 } ;
		    if (isham) list.store (figd) ;
		    OpSum::loop_dt (figd, list) ; 
		    }
		for (int k(0) ;  k < theory.dim ; ++k)
		    {
		    if (i == k || j == k) continue ;
		    Op figa { string {L[i],l[k],l[j],L[k],L[j],l[k],l[i],L[k]}, loop, 4 } ;
		    if (isham) list.store (figa) ;
		    OpSum::loop_dt (figa, list) ; 

		    Op figb { string {L[i],L[k],l[j],l[k],L[j],L[k],l[i],l[k]}, loop, 4 } ;
		    if (isham) list.store (figb) ;
		    OpSum::loop_dt (figb, list) ; 

		    Op figc { string {L[i],l[k],L[j],L[k],l[j],l[k],l[i],L[k]}, loop, 4 } ;
		    if (isham) list.store (figc) ;
		    OpSum::loop_dt (figc, list) ; 

		    Op figd { string {L[i],L[k],L[j],l[k],l[j],L[k],l[i],l[k]}, loop, 4 } ;
		    if (isham) list.store (figd) ;
		    OpSum::loop_dt (figd, list) ; 
		    }
		}
	    }
	for (int i(0) ; i < theory.dim ; ++i)		// polyakov loop
	    {
	    if (theory.box.comp[i])
		{
		Op polyakov  { string (theory.box.comp[i], l[i]), loop, 2 } ;
		if (isham) list.store (polyakov) ;
		OpSum::loop_dt (polyakov, list) ; 

		Op Polyakov  { string (theory.box.comp[i], L[i]), loop, 2 } ;
		if (isham) list.store (Polyakov) ;
		OpSum::loop_dt (Polyakov, list) ; 
		}
	    }
	for (int i(0) ; i < theory.dim ; ++i)		// plaq.polyakov
	    {
	    if (theory.box.comp[i])
		{
		for (int j(0) ; j < theory.dim ; ++j)
		    {
		    if (i == j) continue ;
		    string poly (theory.box.comp[i], l[i]) ;
		    string plaqa {l[i],l[j],L[i],L[j]} ;
		    Op plaqpolya { plaqa+poly, loop, 4 } ;
		    if (isham) list.store (plaqpolya) ;
		    OpSum::loop_dt (plaqpolya, list) ; 

		    string plaqb {l[i],L[j],L[i],l[j]} ;
		    Op plaqpolyb { plaqb+poly, loop, 4 } ;
		    if (isham) list.store (plaqpolyb) ;
		    OpSum::loop_dt (plaqpolyb, list) ; 

		    string Poly (theory.box.comp[i], L[i]) ;
		    string plaqc {L[i],l[j],l[i],L[j]} ;
		    Op plaqPolyc { plaqc+Poly, loop, 4 } ;
		    if (isham) list.store (plaqPolyc) ;
		    OpSum::loop_dt (plaqPolyc, list) ; 

		    string plaqd {L[i],L[j],l[i],l[j]} ;
		    Op plaqPolyd { plaqd+Poly, loop, 4 } ;
		    if (isham) list.store (plaqPolyd) ;
		    OpSum::loop_dt (plaqPolyd, list) ; 
		    }
		}
	    }
	list.setprimary () ;
	}
    else if (theory.nf > 0)
	{
	OpList&	list { global.info(1).ops } ;
	OpType	ferm { OpType::Fermion } ;

	for (int k(0) ; k < theory.nf ; k += 2)		// Fxf, Gxf
	    {
	    for (int i(0) ; i < theory.dim ; ++i)
		{
		Op Gxf { string {F[k+1],l[i],f[k]}, ferm, 1 } ;
		list.store (Gxf) ;
		Op GXf { string {F[k+1],L[i],f[k]}, ferm, 1 } ;
		list.store (GXf) ;
		if (!isham) continue ;
		Op Fxf { string {F[k],l[i],f[k]}, ferm, 1 } ;
		list.store (Fxf) ;
		Op FXf { string {F[k],L[i],f[k]}, ferm, 1 } ;
		list.store (FXf) ;
		}
	    }
	list.setprimary () ;
	}
    }

void OpList::setprimary ()			// Determine Op primacy
    {
    int	opnum ( size() ) ;
    int maxopord (0) ;

    for (auto& op : *this)
	{
	if (!op.order) continue ;
	if (op.order > maxopord) maxopord = op.order ;
	op.primary = true ;
	}
    for (int ord(0) ; ord <= maxopord ; ++ord)
	{
	for (int i(0) ; i < opnum ; ++i)
	    {
	    auto ord1 { (*this)[i].order } ;
	    if (!ord1) continue ;
	    Op op1 { (*this)[i] } ;			// N.B. non-ref needed
	    for (int j(0) ; j <= i ; ++j)		// N.B. include j == i
		{
		auto ord2 { (*this)[j].order } ;
		if (!ord2 || ord1 + ord2 > ord) continue ;
		Op op2 { (*this)[j] } ;			// N.B. non-ref needed
		Gen ans1 (*this) ;

		for (int k(0) ; k <= j ; ++k)		// N.B. include k == j
		    {
		    auto ord3 { (*this)[k].order } ;
		    if (!ord3 || ord1 + ord2 + ord3 != ord) continue ;
		    Op op3 { (*this)[k] } ;		// N.B. non-ref needed

		    if (ans1.empty())
			Commute::op_commute (1.0, op1, op2, ans1) ;

		    for (auto& term : ans1)		// N.B. don't collect
			{
			const Op tmp { (*this)[term.item] } ;
			Gen ans2 (*this) ;
			Commute::op_commute (1.0, op3, tmp, ans2) ;

			for (auto& term : ans2)		// N.B. don't collect
			    {
			    Op& new2 { (*this)[term.item] } ;
			    if (new2.order == ord)
				new2.primary = false ;
			    else if (new2.order > ord)
				cout << "Warning:: mis-ordered Op: " << new2 << "\n" ;
			    }
			}
		    }
		if (ord1 + ord2 != ord) continue ;
		if (ans1.empty())
		    Commute::op_commute (1.0, op1, op2, ans1) ;
		for (auto& term : ans1)			// N.B. don't collect
		    {
		    Op& new1 { (*this)[term.item] } ;
		    if (new1.order == ord)
			new1.primary = false ;
		    else if (new1.order > ord)
			cout << "Warning: mis-ordered Op: " << new1 << "\n" ;
		    }
		}
	    }
	}
    purge (opnum) ;
    }

ObsList::ObsList (const string s, bool can, bool clsfy)	// Construct ObsList
    : name(s), canonicalize(can), classify (clsfy)
    {
    store (Obs(Str(), ObsType::Loop, 0, 0)) ;
    }

PolyTerm ObsList::is_known (Obs& a) const			// Find in ObsList
    {
    int		sgn  ( canonicalize ? a.canon() : a.findstart() ) ;
    numb	indx { find (a) } ;

    return (indx < MAXNUM-1) ? PolyTerm (PolyIndx(indx), sgn)
			     : PolyTerm (PolyIndx(),     0.0) ;
    }

PolyTerm ObsList::is_known (Obs&& a) const			// Find in ObsList
    {
    return is_known (a) ;
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
    const auto&	blab { Blab::level(Blab::OBS) } ;
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
	cout << "catalog " << name << ": stored/found "
	     << a << " -> " << ans[0] << " & "
	     << b << " -> " << ans[1] << "\n" ;
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
	if (size() == MAXNUM-1)
	    gripe ("Max # Obs exceeded: recompile with NUM64!") ;
	push_back (&(iter->first)) ;
	uint stage { o.is_fermi() } ;
	hasher (hash(stage), o) ;
	++nobs (stage) ;
	}
    return indx ;
    }

void ObsList::clear ()					// Empty list
    {
    nobsF = nobsG = 0 ;
    hashF = hashG = 0 ;
    ListMap().swap (map) ;
    vector<const Obs*>().swap (*this) ;
    store (Obs(Str{}, ObsType::Loop, 0, 0)) ;
    }

void ObsList::purge (numb limit)			// Purge entries
    {
    resize (limit) ;
    nobsF -= std::erase_if (map, [&](const auto& p)
	{ return p.second >= limit && p.second >= nobsG; }) ;
    nobsG -= std::erase_if (map, [&](const auto& p)
	{ return p.second >= limit && p.second <  nobsG; }) ;
    if (nobsG + nobsF != size())
	abort ("purge: Inconsistent ObsList size") ;
    rehash () ;
    }

void ObsList::rehash ()				// Recalculate list hashes
    {
    hashF = hashG = 0 ;
    for (const auto& ptr : *this)
	{
	hasher (hash(ptr->is_fermi()), *ptr) ;
	}
    }

void ObsList::hasher (ulong& hash, const Obs& o) const	// List hasher
    {
    ulong x { std::hash<string>{}(o) } ;
    hash ^= x + 0x9e3779b9 + (hash << 6) + (hash >> 2) ;
    }

numb ObsList::find (const Str& s) const			// Find Obs, return index
    {
    auto	p1	{ map.find(s) } ;
    auto&	inbox	{ global.obs.inbox } ;
    if (p1 != map.end()) return p1->second ;
    if (this == &global.obs && inbox.size())
	{
	auto p2 { inbox.find(s) } ;
	if (p2 != inbox.end()) return MAXNUM-1 ;
	}
    return MAXNUM ;
    }

bool ObsList::frozen () const				// Frozen master list?
    {
    return this == &global.obs && global.obs.freeze ;
    }

bool ObsList::freezeif () const				// Freeze if master list
    {
    if (this == &global.obs)
	{
	bool prev { global.obs.freeze } ;
	global.obs.freeze = true ;
	return prev ;
	}
    return false ;
    }

void ObsList::refreezeif (bool prev) const		// Reset freeze if master list
    {
    if (this == &global.obs)
	{
	global.obs.freeze = prev ;
	}
    }

void ObsList::retain (const Obs& o)			// Retain for later insertion
    {
    global.obs.inbox.insert (o) ;
    }


void ObsList::ondisk ()			// Leave global.obs on disk
    {
    auto stage	{ global.stage } ;

    if (size() > 1) cout << "Swapping " << name << " Obs list to disk\n" ;
    if (nobsG > 1)
	{
	bool		is_open	{ global.info(0).sysfile.stream.is_open() } ;
	const auto&	obs	{ global.data(0).obs } ;
	const auto&	obsnum	{ obs.entry().items() } ;
	global.stage = Global::Gauge ;
	if (!is_open || nobsG != obsnum) Save::save_sys () ;
	global.stage = stage ;
	}
    if (nobsF > 0)
	{
	bool		is_open	{ global.info(1).sysfile.stream.is_open() } ;
	const auto&	obs	{ global.data(1).obs } ;
	const auto&	obsnum	{ obs.entry().items() } ;
	global.stage = Global::Fermi ;
	if (!is_open || nobsF != obsnum) Save::save_sys () ;
	global.stage = stage ;
	}
    clear() ;
    global.obs.swapped = true ;
    }

void ObsList::obsinit (uint stage)		// Load basic Obs
    {
    bool iseuc	 { theory.euclid } ;
    char link[4] { 'x', 'y', 'z', 'w' } ;
    char Link[4] { 'X', 'Y', 'Z', 'W' } ;

    global.obs.freeze = false ;
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
	if (!neq (global.obs)) global.info(0).maxord = 2 ;
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
	if (!neq (global.obs))
	    {
	    global.info(1).maxord = 1 ;
	    do_fermiinit() ;
	    }
	}
    global.obs.freeze = true ;

    if (nobsG + nobsF != size())
	abort (format("Inconsistent ObsList size: {} + {} != {}",
	    nobsG, nobsF, size())) ;
    }

int ObsList::do_fermiinit ()			// Initialize fermion -> loop map
    {
    const auto&	blab	{ Blab::level(Blab::OBS) } ;
    uint	initfail (0) ;
    if (blab > 3) cout << "do_fermiinit start\n" << flush ;

    global.obs.fermiinit.clear () ;
    for (numb i(nobsG) ; i < size() ; ++i)
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
	else global.obs.fermiinit.emplace_back (i, term[0]) ;
	}
    return initfail ;
    }

ostream& ObsList::print (ostream& stream, numb indx) const	// Print indexed Obs
    {
    bool	addvev	 { !neq(global.obs) } ;
    auto	prevprec { stream.precision(12) } ;
    const Obs&	obs	 { *at(indx) } ;
    stream << " #" << indx << ": " << std::setw(8) << obs.print() ;
    if (addvev && indx < numerics.vev.size())
	stream << " = " << numerics.vev[indx] ;
    stream << "\n" << setprecision (prevprec) ;
    return stream ;
    }

ostream& ObsList::print (ostream& stream) const			// Print ObsList -> stream
    {
    bool	addvev	 { !neq(global.obs) } ;
    auto	prevprec { stream.precision(12) } ;
    for (int indx(0) ; indx < size() ; ++indx)
	{
	stream << " #" << indx << ": " << *at(indx) ;
	if (addvev && indx < numerics.vev.size())
	    stream << " = " << numerics.vev[indx] ;
	stream << "\n" ;
	}
    stream << setprecision (prevprec) ;
    return stream ;
    }

