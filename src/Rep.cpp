#include "Rep.h"
#include "Gripe.h"
#include <numeric>

Proj::Proj (const string& s, doub d, SymmSum&& v)	// Construct irrep projector
    :
    vector<doub>(Symm::list.size(), 0.0),
    name  (s),
    denom (d)
    {
    for (auto& t : v) (*this)[t.item] += t.coeff ;
    }

Proj Proj::operator() (const Proj& proj) const		// Compose projectors
    {
    Proj	result ;
    for (int i(0) ; i < size() ; ++i)
	{
	doub		coefa = (*this)[i] ;
	const auto&	symma = Symm::list[i] ;
	if (!coefa) continue ;

	for (int j(0) ; j < size() ; ++j)
	    {
	    doub	coefb = proj[j] ;
	    const auto&	symmb = Symm::list[j] ;
	    if (coefb)	result[symma(symmb)] += coefa * coefb ;
	    }
	}
    result.denom = denom * proj.denom ;
    return result ;
    }

Proj Proj::operator+ (const Proj& proj) const		// Add projectors
    {
    Proj sum ;
    int d1 { static_cast<int>(denom) } ;
    int d2 { static_cast<int>(proj.denom) } ;
    int d3 { static_cast<int>(sum.denom) } ;
    if (denom == d1 && proj.denom == d2)
	{
	sum.denom = d3 = std::lcm(d1, d2) ;
	int k1 { d3 / d1 } ;
	int k2 { d3 / d2 } ;

	for (int i(0) ; i < size() ; ++i)
	    {
	    sum[i] = k1 * (*this)[i] + k2 * proj[i] ;
	    }
	}
    else 
	{
	doub k2 = denom / proj.denom ;
	for (int i(0) ; i < size() ; ++i)
	    {
	    sum[i] = (*this)[i] + k2 * proj[i] ;
	    }
	}
    return sum ;
    }

Proj Proj::operator- (const Proj& proj) const		// Subtract projectors
    {
    Proj diff ;
    int d1 { static_cast<int>(denom) } ;
    int d2 { static_cast<int>(proj.denom) } ;
    int d3 { static_cast<int>(diff.denom) } ;
    if (denom == (int) denom && proj.denom == (int) proj.denom)
	{
	diff.denom = d3 = std::lcm(d1, d2) ;
	int k1 { d3 / d1 } ;
	int k2 { d3 / d2 } ;

	for (int i(0) ; i < size() ; ++i)
	    {
	    diff[i] = k1 * (*this)[i] - k2 * proj[i] ;
	    }
	}
    else 
	{
	doub k2 = denom / proj.denom ;
	for (int i(0) ; i < size() ; ++i)
	    {
	    diff[i] = (*this)[i] - k2 * proj[i] ;
	    }
	}
    return diff ;
    }

Proj& Proj::operator+=(const Proj& proj)		// Add to projector
    {
    int d1 { static_cast<int>(denom) } ;
    int d2 { static_cast<int>(proj.denom) } ;
    if (denom == d1 && proj.denom == d2)
	{
	int k1 ( size() / d1 ) ;
	int k2 ( size() / d2 ) ;
	if (k1 != 1)
	    {
	    denom = size() ;
	    for (auto& c : *this) c *= k1 ;
	    }
	for (int i(0) ; i < size() ; ++i)
	    {
	    (*this)[i] += k2 * proj[i] ;
	    }
	}
    else
	{
	doub k2 = denom / proj.denom ;
	for (int i(0) ; i < size() ; ++i)
	    {
	    (*this)[i] += k2 * proj[i] ;
	    }
	}
    return *this ;
    }

Proj& Proj::operator-=(const Proj& proj)		// Subtract from projector
    {
    int d1 { static_cast<int>(denom) } ;
    int d2 { static_cast<int>(proj.denom) } ;
    if (denom == d1 && proj.denom == d2)
	{
	int k1 ( size() / d1 ) ;
	int k2 ( size() / d2 ) ;
	if (k1 != 1)
	    {
	    denom = size() ;
	    for (auto& c : *this) c *= k1 ;
	    }
	for (int i(0) ; i < size() ; ++i)
	    {
	    (*this)[i] -= k2 * proj[i] ;
	    }
	return *this ;
	}
    else
	{
	doub k2 = denom / proj.denom ;
	for (int i(0) ; i < size() ; ++i)
	    {
	    (*this)[i] -= k2 * proj[i] ;
	    }
	}
    return *this ;
    }

bool Proj::allzero () const				// Test if all (essentially) zero
    {
    auto eps { 100 * std::numeric_limits<double>::epsilon() } ;
    return std::all_of (begin(), end(), [eps](doub x) { return abs(x) <= eps ; }) ;
    }

bool Proj::C_even () const				// Charge conjugation even rep?
    {
    for (int i(0) ; i < size()/2 ; ++i)
	{
	auto x { (*this)[i] } ;
	if (!x) continue ;
	auto y { (*this)[i+size()/2] } ;
	if      (y ==  x) return true ;
	else if (y == -x) return false ;
	else fatal ("C indefinite projector") ;
	}
    fatal ("Vanishing projector!") ;
    }

uint2 Proj::indices () const				// Return projector indices
    {
    if (name.size() < 3) return uint2 (0,0) ;

    static string	suffix1 { "abcd" } ;
    static string	suffix2 { "xyzw" } ;
    auto		p { name.end() } ;
    char		j { *--p } ;
    char		i { *--p } ;

    if      (suffix1.find (i) != suffix1.npos)	i -= 'a' ;
    else if (suffix2.find (i) != suffix2.npos)	i -= 'x' ;
    else					i  = 0 ;

    if      (suffix1.find (j) != suffix1.npos)	j -= 'a' ;
    else if (suffix2.find (j) != suffix2.npos)	j -= 'x' ;
    else					j  = 0 ;

    return uint2 (i,j) ;
    }

strview Proj::rowname () const				// Rep row name
    {
    static string suffix { "abcdxyzw" } ;
    if (suffix.find (name.back()) == suffix.npos) return name ;
    else return std::string_view (name.begin(), name.end()-1) ;
    }

strview Proj::repname () const				// Rep name
    {
    static string suffix { "abcdxyzw" } ;
    if (suffix.find (name.back()) == suffix.npos) return name ;
    else return std::string_view (name.begin(), name.end()-2) ;
    }

ostream& operator<< (ostream& stream, const Proj& proj)	// Print projector
    {
    stream << proj.name << " = 1/" << proj.denom << " (" ;
    for (int i(0) ; i < proj.size() ; ++i)
	{
	doub coeff { proj[i] } ;
	if (coeff != 0)
	    {
	    string suffix ;
	    if      (coeff ==  1.0)	stream << " +" ;
	    else if (coeff == -1.0)	stream << " -" ;
	    else if (coeff ==  0.5)	stream << " +", suffix = "/2" ;
	    else if (coeff == -0.5)	stream << " -", suffix = "/2" ;
	    else			stream << format(" {:+.1f} ", coeff) ;
	    stream << Symm::list[i].name << suffix ;
	    }
	}
    stream << " )\n" ;
    return stream ;
    }

uint Rep::known (const string& name)			// Return named rep index
    {
    return Rep::list.map.at(name) ;
    }

ostream& operator<< (ostream& stream, const Rep& rep)	// Print Rep
    {
    for (auto& proj : rep) stream << proj ;
    return stream ;
    }

void Rep::repinit()					// Define symmetry projectors
    {
    double	sq3  { sqrt(3) } ;
    auto&	proj { Proj::list } ;
    SymmTerm	(*s)(const string&&) = Symm::known ;

    if constexpr (theory.isR1() || theory.isS1())	// 1D lattice
	{
	proj.emplace_back("A1p", 4,     (s("E") + s("C")) *
					(s("E") + s("Rx")) ) ;

	proj.emplace_back("A2p", 4,     (s("E") + s("C")) *
					(s("E") - s("Rx")) ) ;

	proj.emplace_back("A1m", 4,     (s("E") - s("C")) *
					(s("E") + s("Rx")) ) ;

	proj.emplace_back("A2m", 4,     (s("E") - s("C")) *
					(s("E") - s("Rx")) ) ;
	}
    else if constexpr (theory.isR2())			// 2D lattice
	{
	proj.emplace_back("A1p", 16,    (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("A2p", 16,    (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B1p", 16,    (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B2p", 16,    (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("Epxx",  8,   (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) ) ;

	proj.emplace_back("Epxy",  8,   (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Epyx",  8,   (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Epyy",  8,   (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) ) ;

	proj.emplace_back("A1m", 16,    (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("A2m", 16,    (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B1m", 16,    (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B2m", 16,    (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("Emxx",  8,   (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) ) ;

	proj.emplace_back("Emxy",  8,   (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Emyx",  8,   (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Emyy",  8,   (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) ) ;

	}
    else if constexpr (theory.isR3())			// 3D lattice
	{
	proj.emplace_back("A1gp", 96,   (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxyz") + s("Pxzy")
					+ s("Pxy") + s("Pxz") + s("Pyz")) ) ;

	proj.emplace_back("A2gp", 96,   (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxyz") + s("Pxzy")
					- s("Pxy") - s("Pxz") - s("Pyz")) ) ;

	proj.emplace_back("Egpaa", 48,  (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxyz")/2 - s("Pxzy")/2
					+ s("Pxy") - s("Pxz")/2 - s("Pyz")/2) ) ;

	proj.emplace_back("Egpab", 32*sq3, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("Pxyz") - s("Pxzy") - s("Pxz") + s("Pyz")) ) ;

	proj.emplace_back("Egpba", 32*sq3, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("Pxzy") - s("Pxyz") - s("Pxz") + s("Pyz")) ) ;

	proj.emplace_back("Egpbb", 48,  (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxyz")/2 - s("Pxzy")/2
					- s("Pxy") + s("Pxz")/2 + s("Pyz")/2) ) ;

	proj.emplace_back("T1gpxx", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pyz")) ) ;

	proj.emplace_back("T1gpxy", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pyz")) * s("Pxy") ) ;

	proj.emplace_back("T1gpxz", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pyz")) * s("Pxz") ) ;

	proj.emplace_back("T1gpyx", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxz")) * s("Pxy") ) ;

	proj.emplace_back("T1gpyy", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxz")) ) ;

	proj.emplace_back("T1gpyz", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("Pxz") - s("E")) * s("Pyz") ) ;

	proj.emplace_back("T1gpzx", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxy")) * s("Pxz") ) ;

	proj.emplace_back("T1gpzy", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("Pxy") - s("E")) * s("Pyz") ) ;

	proj.emplace_back("T1gpzz", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("T2gpaa", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pyz")) ) ;

	proj.emplace_back("T2gpab", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pyz")) * s("Pxy") ) ;

	proj.emplace_back("T2gpac", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pyz")) * s("Pxz") ) ;

	proj.emplace_back("T2gpba", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxz")) * s("Pxy") ) ;

	proj.emplace_back("T2gpbb", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxz")) ) ;

	proj.emplace_back("T2gpbc", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxz")) * s("Pyz") ) ;

	proj.emplace_back("T2gpca", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxy")) * s("Pxz") ) ;

	proj.emplace_back("T2gpcb", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxy")) * s("Pyz") ) ;

	proj.emplace_back("T2gpcc", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("A1up", 96,   (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxyz") + s("Pxzy")
					- s("Pxy") - s("Pxz") - s("Pyz")) ) ;

	proj.emplace_back("A2up", 96,   (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxyz") + s("Pxzy")
					+ s("Pxy") + s("Pxz") + s("Pyz")) ) ;

	proj.emplace_back("Eupaa", 48,  (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxyz")/2 - s("Pxzy")/2
					- s("Pxy") + s("Pxz")/2 + s("Pyz")/2) ) ;

	proj.emplace_back("Eupab", 32*sq3, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("Pxyz") - s("Pxzy") + s("Pxz") - s("Pyz")) ) ;

	proj.emplace_back("Eupba", 32*sq3, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("Pxzy") - s("Pxyz") + s("Pxz") - s("Pyz")) ) ;

	proj.emplace_back("Eupbb", 48,  (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxyz")/2 - s("Pxzy")/2
					+ s("Pxy") - s("Pxz")/2 - s("Pyz")/2) ) ;

	proj.emplace_back("T1upxx", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pyz")) ) ;

	proj.emplace_back("T1upxy", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pyz")) * s("Pxy") ) ;

	proj.emplace_back("T1upxz", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pyz")) * s("Pxz") ) ;

	proj.emplace_back("T1upyx", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxz")) * s("Pxy") ) ;

	proj.emplace_back("T1upyy", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxz")) ) ;

	proj.emplace_back("T1upyz", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxz")) * s("Pyz") ) ;

	proj.emplace_back("T1upzx", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxy")) * s("Pxz") ) ;

	proj.emplace_back("T1upzy", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxy")) * s("Pyz") ) ;

	proj.emplace_back("T1upzz", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("T2upaa", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pyz")) ) ;

	proj.emplace_back("T2upab", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pyz")) * s("Pxy") ) ;

	proj.emplace_back("T2upac", 32, (s("E") + s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pyz")) * s("Pxz") ) ;

	proj.emplace_back("T2upba", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxz")) * s("Pxy") ) ;

	proj.emplace_back("T2upbb", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxz")) ) ;

	proj.emplace_back("T2upbc", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("Pxz") - s("E")) * s("Pyz") ) ;

	proj.emplace_back("T2upca", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxy")) * s("Pxz") ) ;

	proj.emplace_back("T2upcb", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("Pxy") - s("E")) * s("Pyz") ) ;

	proj.emplace_back("T2upcc", 32, (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("A1gm", 96,   (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxyz") + s("Pxzy")
					+ s("Pxy") + s("Pxz") + s("Pyz")) ) ;

	proj.emplace_back("A2gm", 96,   (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxyz") + s("Pxzy")
					- s("Pxy") - s("Pxz") - s("Pyz")) ) ;

	proj.emplace_back("Egmaa", 48,  (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxyz")/2 - s("Pxzy")/2
					+ s("Pxy") - s("Pxz")/2 - s("Pyz")/2) ) ;

	proj.emplace_back("Egmab", 32*sq3, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("Pxyz") - s("Pxzy") - s("Pxz") + s("Pyz")) ) ;

	proj.emplace_back("Egmba", 32*sq3, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("Pxzy") - s("Pxyz") - s("Pxz") + s("Pyz")) ) ;

	proj.emplace_back("Egmbb", 48,  (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxyz")/2 - s("Pxzy")/2
					- s("Pxy") + s("Pxz")/2 + s("Pyz")/2) ) ;

	proj.emplace_back("T1gmxx", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pyz")) ) ;

	proj.emplace_back("T1gmxy", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pyz")) * s("Pxy") ) ;

	proj.emplace_back("T1gmxz", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pyz")) * s("Pxz") ) ;

	proj.emplace_back("T1gmyx", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxz")) * s("Pxy") ) ;

	proj.emplace_back("T1gmyy", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxz")) ) ;

	proj.emplace_back("T1gmyz", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("Pxz") - s("E")) * s("Pyz") ) ;

	proj.emplace_back("T1gmzx", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxy")) * s("Pxz") ) ;

	proj.emplace_back("T1gmzy", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("Pxy") - s("E")) * s("Pyz") ) ;

	proj.emplace_back("T1gmzz", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("T2gmaa", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pyz")) ) ;

	proj.emplace_back("T2gmab", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pyz")) * s("Pxy") ) ;

	proj.emplace_back("T2gmac", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pyz")) * s("Pxz") ) ;

	proj.emplace_back("T2gmba", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxz")) * s("Pxy") ) ;

	proj.emplace_back("T2gmbb", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxz")) ) ;

	proj.emplace_back("T2gmbc", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxz")) * s("Pyz") ) ;

	proj.emplace_back("T2gmca", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxy")) * s("Pxz") ) ;

	proj.emplace_back("T2gmcb", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxy")) * s("Pyz") ) ;

	proj.emplace_back("T2gmcc", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("A1um", 96,   (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxyz") + s("Pxzy")
					- s("Pxy") - s("Pxz") - s("Pyz")) ) ;

	proj.emplace_back("A2um", 96,   (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxyz") + s("Pxzy")
					+ s("Pxy") + s("Pxz") + s("Pyz")) ) ;

	proj.emplace_back("Eumaa", 48,  (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxyz")/2 - s("Pxzy")/2
					- s("Pxy") + s("Pxz")/2 + s("Pyz")/2) ) ;

	proj.emplace_back("Eumab", 32*sq3, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("Pxyz") - s("Pxzy") + s("Pxz") - s("Pyz")) ) ;

	proj.emplace_back("Eumba", 32*sq3, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("Pxzy") - s("Pxyz") + s("Pxz") - s("Pyz")) ) ;

	proj.emplace_back("Eumbb", 48,  (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxyz")/2 - s("Pxzy")/2
					+ s("Pxy") - s("Pxz")/2 - s("Pyz")/2) ) ;

	proj.emplace_back("T1umxx", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pyz")) ) ;

	proj.emplace_back("T1umxy", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pyz")) * s("Pxy") ) ;

	proj.emplace_back("T1umxz", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pyz")) * s("Pxz") ) ;

	proj.emplace_back("T1umyx", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxz")) * s("Pxy") ) ;

	proj.emplace_back("T1umyy", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxz")) ) ;

	proj.emplace_back("T1umyz", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Pxz")) * s("Pyz") ) ;

	proj.emplace_back("T1umzx", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxy")) * s("Pxz") ) ;

	proj.emplace_back("T1umzy", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxy")) * s("Pyz") ) ;

	proj.emplace_back("T1umzz", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("T2umaa", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pyz")) ) ;

	proj.emplace_back("T2umab", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pyz")) * s("Pxy") ) ;

	proj.emplace_back("T2umac", 32, (s("E") - s("C")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pyz")) * s("Pxz") ) ;

	proj.emplace_back("T2umba", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxz")) * s("Pxy") ) ;

	proj.emplace_back("T2umbb", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Pxz")) ) ;

	proj.emplace_back("T2umbc", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Rz")) *
					(s("Pxz") - s("E")) * s("Pyz") ) ;

	proj.emplace_back("T2umca", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxy")) * s("Pxz") ) ;

	proj.emplace_back("T2umcb", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("Pxy") - s("E")) * s("Pyz") ) ;

	proj.emplace_back("T2umcc", 32, (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Pxy")) ) ;

	}
    else if constexpr (theory.isR4())			// 4D lattice
	{
	proj.emplace_back("Ap", 768,    (s("E") + s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rw")) *
					(s("E")
					+ s("Pxy") + s("Pxz") + s("Pxw")
					+ s("Pyz") + s("Pyw") + s("Pzw")
					+ s("Pxyz") + s("Pxyw")
					+ s("Pxzy") + s("Pxzw")
					+ s("Pxwy") + s("Pxwz")
					+ s("Pyzw") + s("Pywz")
					+ s("Pxyzw") + s("Pxywz")
					+ s("Pxzyw") + s("Pxzwy")
					+ s("Pxwyz") + s("Pxwzy")
					+ s("PxyPzw") + s("PxzPyw") + s("PxwPyz")) ) ;

	proj.emplace_back("Am", 768,    (s("E") - s("C")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rw")) *
					(s("E")
					+ s("Pxy") + s("Pxz") + s("Pxw")
					+ s("Pyz") + s("Pyw") + s("Pzw")
					+ s("Pxyz") + s("Pxyw")
					+ s("Pxzy") + s("Pxzw")
					+ s("Pxwy") + s("Pxwz")
					+ s("Pyzw") + s("Pywz")
					+ s("Pxyzw") + s("Pxywz")
					+ s("Pxzyw") + s("Pxzwy")
					+ s("Pxwyz") + s("Pxwzy")
					+ s("PxyPzw") + s("PxzPyw") + s("PxwPyz")) ) ;
	}
    else if constexpr (theory.isR2xS1())	// 3D lattice, S(1) compactified
	{
	proj.emplace_back("A1gp", 32,   (s("E") + s("C")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("A2gp", 32,   (s("E") + s("C")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B1gp", 32,   (s("E") + s("C")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B2gp", 32,   (s("E") + s("C")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("Egpxx", 16,  (s("E") + s("C")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) ) ;

	proj.emplace_back("Egpxy", 16,  (s("E") + s("C")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Egpyx", 16,  (s("E") + s("C")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Egpyy", 16,  (s("E") + s("C")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) ) ;

	proj.emplace_back("A1up", 32,   (s("E") + s("C")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("A2up", 32,   (s("E") + s("C")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B1up", 32,   (s("E") + s("C")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B2up", 32,   (s("E") + s("C")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("Eupxx", 16,  (s("E") + s("C")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) ) ;

	proj.emplace_back("Eupxy", 16,  (s("E") + s("C")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Eupyx", 16,  (s("E") + s("C")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Eupyy", 16,  (s("E") + s("C")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) ) ;

	proj.emplace_back("A1gm", 32,   (s("E") - s("C")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("A2gm", 32,   (s("E") - s("C")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B1gm", 32,   (s("E") - s("C")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B2gm", 32,   (s("E") - s("C")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("Egmxx", 16,  (s("E") - s("C")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) ) ;

	proj.emplace_back("Egmxy", 16,  (s("E") - s("C")) *
					(s("E") + s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Egmyx", 16,  (s("E") - s("C")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Egmyy", 16,  (s("E") - s("C")) *
					(s("E") + s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) ) ;

	proj.emplace_back("A1um", 32,   (s("E") - s("C")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("A2um", 32,   (s("E") - s("C")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B1um", 32,   (s("E") - s("C")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") + s("Ry")) *
					(s("E") - s("Pxy")) ) ;

	proj.emplace_back("B2um", 32,   (s("E") - s("C")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") - s("Ry")) *
					(s("E") + s("Pxy")) ) ;

	proj.emplace_back("Eumxx", 16,  (s("E") - s("C")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) ) ;

	proj.emplace_back("Eumxy", 16,  (s("E") - s("C")) *
					(s("E") - s("Rz")) *
					(s("E") - s("Rx")) *
					(s("E") + s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Eumyx", 16,  (s("E") - s("C")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) * s("Pxy") ) ;

	proj.emplace_back("Eumyy", 16,  (s("E") - s("C")) *
					(s("E") - s("Rz")) *
					(s("E") + s("Rx")) *
					(s("E") - s("Ry")) ) ;

	}
    for (int i(0) ; i < proj.size() ; ++i)
	{
	string	rep { proj[i].rowname() } ;
	auto	[iter, isnew] { list.map.try_emplace (rep, list.size()) } ;

	if (isnew) list.push_back (Rep(rep, proj[i].C_even())) ;
	list[iter->second].push_back (proj[i]) ;
	}
    if (theory.nrep != list.size()) fatal ("theory.nrep is wrong!") ;
    }
