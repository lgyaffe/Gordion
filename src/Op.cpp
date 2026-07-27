#include "Op.h"
#include "Commute.h"
#include "Global.h"
#include "Gripe.h"
#include <numeric>
#include <regex>

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
    else if (fs && !isF(front()))	err = "is malformed fermion bilinear" ;
    else if (fs && !isf(back()))	err = "is malformed fermion bilinear" ;
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

OpSum::OpSum (OpList& oplist)			// Constructor
    :
    vector<OpTerm>::vector(), list (oplist)
    {}

OpSum::OpSum (OpTerm* beg, OpTerm* end, OpList& oplist)
    :
    vector<OpTerm>::vector(beg, end), list (oplist)
    {}

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
    OpSum ans { list } ;;
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

void OpList::setprimary ()			// Determine Op primacy
    {
    int	opnum	( size() ) ;
    int maxord	(0) ;

    for (int i(0) ; i < opnum ; ++i)
	{
	Op& op { (*this)[i] } ;
	if (!op.order) continue ;
	if (op.order > maxord) maxord = op.order ;
	op.primary = true ;
	}
    for (int i(0) ; i < opnum ; ++i)
	{
	Op op1 { (*this)[i] } ;			// N.B. non-ref needed
	if (!op1.order) continue ;
	for (int j(0) ; j <= i ; ++j)		// N.B. include i == j
	    {
	    Op op2 { (*this)[j] } ;		// N.B. non-ref needed
	    if (!op2.order) continue ;
	    if (op1.order + op2.order > maxord)	continue ;
	    Gen ans1 (*this) ;
	    Commute::op_commute (1.0, op1, op2, ans1) ;

	    for (auto& term : ans1)		// N.B. don't collect
		{
		Op& new1 { (*this)[term.item] } ;
		if (new1.order == op1.order + op2.order)
		    new1.primary = false ;
		else if (new1.order > op1.order + op2.order)
		    cout << "Warning: mis-ordered Op: " << new1 << "\n" ;
		}
	    for (int k(0) ; k <= i ; ++k)	// need triples for fermions
		{
		Op op3 { (*this)[k] } ;		// N.B. non-ref needed
		if (!op3.order) continue ;
		if (op1.order + op2.order + op3.order > maxord) continue ;
		for (auto& term : ans1)		// N.B. don't collect
		    {
		    const Op tmp { (*this)[term.item] } ;
		    Gen ans2 (*this) ;
		    Commute::op_commute (1.0, op3, tmp, ans2) ;

		    for (auto& term : ans2)	// N.B. don't collect
			{
			Op& new2 { (*this)[term.item] } ;
			if (new2.order == op1.order + op2.order + op3.order)
			    new2.primary = false ;
			else if (new2.order > op1.order + op2.order + op3.order)
			    cout << "Warning: mis-ordered Op: " << new2 << "\n" ;
			}
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

