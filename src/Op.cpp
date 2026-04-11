#include "Op.h"
#include "Global.h"
#include "Commute.h"
#include "Gripe.h"
#include <numeric>
#include <regex>

Op::Op(const string& s, OpType t, short ord)	// Construct from string
    :
    SymbStr(to_symb(s)), type(t), order(ord)
    {
    if (type == OpType::Loop) findstart() ;
    validate() ;
    }

Op::Op(const string& s, short ord)		// Construct from string
    :
    SymbStr(to_symb(s)), type(optype(s)), order(ord)
    {
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

    if (err) gripe (format("Bad Op: {} {}", SymbStr::print(), err)) ;
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

void Op::setprimary ()				// Determine Op primacy
    {
    int  stage	{ global.stage } ;
    auto maxgen { global.info().maxgen } ;
    int	 opnum	( list.size() ) ;

    for (int i(0) ; i < opnum ; ++i)
	{
	const Op& op { list[i] } ;
	if (op.is_Fermion() == stage && op.order <= maxgen) op.primary = true ;
	}
    for (int i(0) ; i < opnum ; ++i)
	{
	Op op1 { list[i] } ;		// N.B. non-ref 'cuz list may grow
	if (!op1.order || op1.is_Fermion() != stage) continue ;

	for (int j(0) ; j < i ; ++j)
	    {
	    Op op2 { list[j] } ;	// N.B. non-ref 'cuz list may grow
	    if (!op2.order || op2.is_Fermion() != stage) continue ;

	    if (op1.order + op2.order <= global.info().maxgen)
		{
		Gen ans ;
		Commute::op_commute (1.0, op1, op2, ans) ;

		for (auto& term : ans)
		    {
		    const Op& op3 { list[term.item] } ;
		    if (op3.order == op1.order + op2.order)
			{
			op3.primary = false ;
			}
		    }
		}
	    }
	}
    Op::purge (opnum) ;
    }

OpSum Op::loop_dt (Op op)			// Loop Op -> Eloop OpSum
    {
    if (op.type != OpType::Loop) fatal ("Bad call to loop_dt") ;
    OpSum ans ;
    return loop_dt (OpTerm(Op::store(op)), ans) ;
    }

OpSum Op::loop_dt (OpSum s)			// Loop OpSum -> Eloop OpSum
    {
    OpSum ans ;
    for (auto& t : s) loop_dt (t, ans) ;
    return ans ;
    }

OpSum Op::flipT (OpSum s)			// Flip bilinear staggering
    {
    cout << "Op::flipT OpSum: " << s.size() << "\n" ;
    OpSum ans ;
    for (auto& t : s)
	{
	cout << "t " << list[t.item] << "\n" ;
	Op op { list[t.item] } ;
	if (op.type != OpType::Fermion) fatal ("Bad call to flipT") ;
	op.front() = stag(op.front()) ;
	cout << " storing " << op << "\n" ;
	ans.emplace_back ( Op::store(op) ) ;
	cout << " ans.size " << ans.size() << "\n" ;
	}
    return ans ;
    }

OpSum Op::loop_dt (OpTerm t, OpSum& ans)	// Loop OpTerm -> Eloop OpSum
    {
    Op op { list[t.item] } ;
    if (op.type != OpType::Loop) fatal ("Bad call to loop_dt") ;
    op.type = OpType::Eloop ;

    for (auto ptr = op.begin() ; ptr < op.end() ; ++ptr)
	{
	op.front() += addE ;
	uint	indx { Op::store (op) } ;
	doub	coef { isrefl(op.front()) ? -t.coeff : t.coeff } ;
	ans.emplace_back ( Op::store(op), coef ) ;
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
    resize(distance(begin(),a)) ;

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

uint Op::store (const Op& op)				// Store Op in list
    {
    uint	nop  ( list.size() ) ;
    uint	indx { list.store (op) } ;

    if (list.size() > nop)		// Op added
	{
	++global.info (op.is_Fermion()).nop ;
	if (global.nopF() && !op.is_Fermion())
	    gripe ("Can't add gauge Op after fermion Op's") ;
	}
    if (indx >= list.size()) fatal ("Op::store: bad store! ") ;
    return indx ;
    }

ostream& Op::print (ostream& stream, uint indx)		// Print indexed Op
    {
    const Op& op { list[indx] } ;
    return stream << " op #" << indx << " = " << op << "\n" ;
    }

ostream& Op::print (ostream& stream)			// Print Op::list
    {
    stream << " operators:\n" ;
    for (int indx(0) ; indx < list.size() ; ++indx)
	{
	const Op& op { list[indx] } ;
	stream  << " #" << indx << " = " << op << "\n" ;
	}
    return stream ;
    }

ostream& operator<< (ostream& stream, const Op& op)	// Print Op
    {
    if (op.size()) stream << static_cast<SymbStr>(op) ;
    else stream << "1" ;
    if (op.order >= 0) stream << " (" << op.order << ")" ;
    return stream ;
    }

ostream& operator<< (ostream& stream, const OpSum& s)	// Print OpSum
    {
    for (auto& t : s)
	{
	coeffprt (stream, t.coeff) ;
	stream << Op::list[t.item] ;
	}
    return stream ;
    }
