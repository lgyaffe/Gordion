#include "Op.h"
#include "Commute.h"
#include "Global.h"
#include "Gripe.h"
#include <numeric>
#include <regex>

OpSum::OpSum (OpList& oplist)			// Constructor
    :
    vector<OpTerm>::vector(), list (oplist)
    {}

OpSum::OpSum (OpTerm* beg, OpTerm* end, OpList& oplist)
    :
    vector<OpTerm>::vector(beg, end), list (oplist)
    {}

Op::Op(const string& s, OpType t, short ord)	// Construct from string
    :
    Str(s), type(t), order(ord)
    {
    if (order > MAXORD) gripe ("Max Op order exceeded: recompile without NUM32") ;
    if (type == OpType::Loop) findstart() ;
    joinends() ;
    validate() ;
    }

Op::Op(const string& s, short ord)		// Construct from string
    :
    Str(s), type(optype(s)), order(ord)
    {
    if (order > MAXORD) gripe ("Max Op order exceeded: recompile without NUM32") ;
    joinends() ;
    if (type == OpType::Loop) findstart() ;
    validate() ;
    }

void Op::validate()				// Test Op validity
    {
    int  fs	(0) ;
    int  Fs	(0) ;
    int  Es	(0) ;
    int  nf	(0) ;
    int  dim	(0) ;
    int  derivs	(0) ;
    bool iseuc	{ theory.euclid } ;

    const char* err = nullptr ;
    for (auto c : *this)
	{
	if (islink(c) && dim <= axis(c))	dim = axis(c)+1 ;
	if (isferm(c) && nf  <= flav(c))	nf  = flav(c)/2+1 ;
	if (isferm(c) && !isconj(c))		++fs ;
	if (isferm(c) &&  isconj(c))		++Fs ;
	if (isferm(c) && isderiv(c))		++derivs ;
	if (isE(c)  || isElink(c))		++derivs, ++Es ;
	if (isEE(c) || isEElink(c))		Es += 2 ;
	}
    if      (nf  > theory.nf)		err = "has excess fermion flavors" ;
    else if (dim > theory.dim)		err = "exceeds lattice dimension" ;
    else if ( fs && (fs != Fs))		err = "has bad fermion insertions" ;
    else if (!fs && !isclosed())	err = "is not closed loop" ;
    else if (fs && !is_F(front()))	err = "is malformed fermion bilinear" ;
    else if (fs && !is_f(back()))	err = "is malformed fermion bilinear" ;
    else if (fs > 1)			err = "has excessive fermions" ;
    else if (Es && islink(front()))	err = "is mis-rotated" ;
    else if (Es + fs > 1)		err = "has excessive E's" ;
    else if (iseuc && derivs > 1)	err = "has too many derivatives" ;
    else if (iseuc && Fs && !derivs)	err = "has no derivative" ;
    else if (fs && type != OpType::Fermion)	  err = "has wrong type" ;
    else if (!fs && !Es && type != OpType::Loop)  err = "has wrong type" ;
    else if (!fs &&  Es && type != OpType::Eloop) err = "has wrong type" ;

    if (err) gripe (format("Bad Op: {} {}", Str::print(), err)) ;
    }

OpType Op::optype (const string s)		// Determine Op type
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

    if      (nE == 0 && nF == 0) return OpType::Loop ;
    else if (nE == 1 && nF == 0) return OpType::Eloop ;
    else if (nE == 0 && nF == 1) return OpType::Fermion ;
    else gripe (format ("Op {} is invalid type", s)) ;
    }

void Op::findstart()				// Rotate to preferred start
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

OpSum OpSum::flipT () const		// Flip bilinear staggering
    {
    OpSum ans { oplist() } ;
    for (auto& t : *this)
	{
	Op op { oplist()[t.item] } ;
	if (op.type != OpType::Fermion) fatal ("Bad call to flipT") ;
	op.front() = stag(op.front()) ;
	ans.emplace_back ( ans.oplist().store(op) ) ;
	}
    return ans ;
    }

OpSum OpSum::loop_dt ()				// Loop OpSum -> Eloop OpSum
    {
    OpSum ans { oplist() } ;
    for (auto& t : *this) loop_dt (t, ans) ;
    return ans ;
    }

OpSum OpSum::loop_dt (Op op, OpList& list)	// Loop Op -> Eloop OpSum
    {
    if (op.type != OpType::Loop) fatal ("Bad call to loop_dt") ;
    OpSum ans { list } ;
    return loop_dt (OpTerm (list.store(op)), ans) ;
    }

OpSum OpSum::loop_dt (OpTerm t, OpSum& ans)	// Loop OpTerm -> Eloop OpSum
    {
    Op op { ans.oplist()[t.item] } ;
    if (op.type != OpType::Loop) fatal ("Bad call to loop_dt") ;
    op.type = OpType::Eloop ;

    for (auto ptr = op.begin() ; ptr < op.end() ; ++ptr)
	{
	op.front() += addE ;
	numb	indx { ans.oplist().store (op) } ;
	real	coef { isrefl(op.front()) ? -t.coeff : t.coeff } ;
	ans.emplace_back ( indx, coef ) ;
	op.front() -= addE ;
	rotate (op.begin(), op.begin() + 1, op.end()) ;
	}
    return ans ;
    }

int OpSum::collect (bool divgcd)		// Collect terms, optionally
    {						// divide by & return gcd
    std::sort(begin(), end(),
	[](const OpTerm& a, const OpTerm& b) { return a.item < b.item ; });

    auto a = begin() ;
    for (auto b = begin() ; b < end() ; ++a)
	{
	if (b > a) *a = *b ;
	while (++b < end() && b->item == a->item) a->coeff += b->coeff ;
	}
    resize (distance(begin(),a)) ;

    if (divgcd)
	{
	int k(0) ;
	for (auto a = begin() ; a < end() ; ++a)
	    {
	    int j = static_cast<int>(a->coeff) ;
	    if (a->coeff == j) k = std::gcd(j,k) ;
	    else return 1 ;
	    }
	if (k > 1) for (auto a = begin() ; a < end() ; ++a) a->coeff /= k ;
	return k ;
	}
    else return 1 ;
    }

ostream& operator<< (ostream& stream, const Op& op)	// Print Op
    {
    if (op.size()) stream << static_cast<Str>(op) ;
    else stream << "1" ;
    if (op.order >= 0) stream << " (" << op.order << ")" ;
    return stream ;
    }

ostream& operator<< (ostream& stream, const OpSum& s)	// Print OpSum
    {
    for (auto& t : s)
	{
	Print::coeffprt (stream, t.coeff) ;
	stream << s.oplist() [t.item] ;
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
	    for (int j(i) ; ++j < theory.dim ;)
		{
		Op plaq { string {l[i],l[j],L[i],L[j]}, loop, 2 } ;
		Op Plaq { string {l[i],L[j],L[i],l[j]}, loop, 2 } ;
		if (isham) list.store (plaq) ;
		if (isham) list.store (Plaq) ;
		OpSum::loop_dt (plaq, list) ;
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
		    Op fig8 { string {l[i],l[k],l[j],L[k],L[j],l[k],L[i],L[k]}, loop, 4 } ;
		    Op Fig8 { string {l[i],L[k],l[j],l[k],L[j],L[k],L[i],l[k]}, loop, 4 } ;
		    if (isham) list.store (fig8) ;
		    if (isham) list.store (Fig8) ;
		    OpSum::loop_dt (fig8, list) ; 
		    OpSum::loop_dt (Fig8, list) ; 
		    }
		}
	    for (int j(0) ;  j < theory.dim ; ++j)
		{
		for (int k(0) ;  k < theory.dim ; ++k)
		    {
		    if (i == k || j == k) continue ;
		    Op fig8 { string {l[i],l[k],L[j],L[k],l[j],l[k],L[i],L[k]}, loop, 4 } ;
		    Op Fig8 { string {l[i],L[k],L[j],l[k],l[j],L[k],L[i],l[k]}, loop, 4 } ;
		    if (isham) list.store (fig8) ;
		    if (isham) list.store (Fig8) ;
		    OpSum::loop_dt (fig8, list) ; 
		    OpSum::loop_dt (Fig8, list) ; 
		    }
		}
	    }
	for (int i(0) ; i < theory.dim ; ++i)		// polyakov loop
	    {
	    if (theory.box.comp[i])
		{
		Op polyakov  { string (theory.box.comp[i], l[i]), loop, 2 } ;
		Op Polyakov  { string (theory.box.comp[i], L[i]), loop, 2 } ;
		if (isham) list.store (polyakov) ;
		if (isham) list.store (Polyakov) ;
		OpSum::loop_dt (polyakov, list) ; 
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
		    string plaq {l[i],l[j],L[i],L[j]} ;
		    string Plaq {l[i],L[j],L[i],l[j]} ;
		    string poly (theory.box.comp[i], l[i]) ;
		    Op plaqpoly { plaq+poly, loop, 4 } ;
		    Op Plaqpoly { Plaq+poly, loop, 4 } ;
		    if (isham) list.store (plaqpoly) ;
		    if (isham) list.store (Plaqpoly) ;
		    OpSum::loop_dt (plaqpoly, list) ; 
		    OpSum::loop_dt (Plaqpoly, list) ; 
		    }
		for (int j(0) ; j < theory.dim ; ++j)
		    {
		    if (i == j) continue ;
		    string plaq {L[i],l[j],l[i],L[j]} ;
		    string Plaq {L[i],L[j],l[i],l[j]} ;
		    string Poly (theory.box.comp[i], L[i]) ;
		    Op plaqPoly { plaq+Poly, loop, 4 } ;
		    Op PlaqPoly { Plaq+Poly, loop, 4 } ;
		    if (isham) list.store (plaqPoly) ;
		    if (isham) list.store (PlaqPoly) ;
		    OpSum::loop_dt (plaqPoly, list) ; 
		    OpSum::loop_dt (PlaqPoly, list) ; 
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
		Op GXf { string {F[k+1],L[i],f[k]}, ferm, 1 } ;
		list.store (Gxf) ;
		list.store (GXf) ;
		if (!isham) continue ;
		Op Fxf { string {F[k],l[i],f[k]}, ferm, 1 } ;
		Op FXf { string {F[k],L[i],f[k]}, ferm, 1 } ;
		list.store (Fxf) ;
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
		    Op op3 { (*this)[j] } ;		// N.B. non-ref needed

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

