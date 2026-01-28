#include "Symm.h"
#include "Blab.h"
#include <sstream>
#include <cstring>

const array<symb,Nsymb> Symm::idmap				// Identity map
    {
    '\x00','\x01','\x02','\x03','\x04','\x05','\x06','\x07',
    '\x08','\x09','\x0a','\x0b','\x0c','\x0d','\x0e','\x0f',
    '\x10','\x11','\x12','\x13','\x14','\x15','\x16','\x17',
    '\x18','\x19','\x1a','\x1b','\x1c','\x1d','\x1e','\x1f',
    '\x20','\x21','\x22','\x23','\x24','\x25','\x26','\x27',
    '\x28','\x29','\x2a','\x2b','\x2c','\x2d','\x2e','\x2f'
    } ;

std::size_t Symmhash::operator()(const Symm& s) const noexcept	// Symm hash function
    {
    string buf(s.map.size()+1, '\0') ;
    std::memcpy (&buf[0], s.map.begin(), s.map.size()) ;
    std::size_t h1 = std::hash<std::string>{}(buf) ;
    std::size_t h2 = std::hash<std::bitset<Nsymb+1>>{}(s.sgn) ;
    return h1 ^ (h2 << 1) ;
    }

SymmSum SymmTerm::operator+(SymmTerm a) const			// Add two symmetries
    {
    return SymmSum { *this, a } ;
    }

SymmSum SymmTerm::operator-(SymmTerm a) const			// Subtract symmetries
    {
    return SymmSum { *this, SymmTerm { a.item, -a.coeff } } ;
    }

SymmSum& SymmSum::operator+(SymmTerm a)				// Add symmetry term to sum
    {
    this->push_back(a) ;
    return *this ;
    } ;

SymmSum& SymmSum::operator-(SymmTerm a)				// Subtract term from sum
    {
    this->emplace_back (a.item, -a.coeff) ;
    return *this ;
    } ;

SymmSum SymmSum::operator*(SymmSum v)				// Multiply symmetry sums
    {
    SymmSum prod ;
    for (const auto& a : *this)
	{
	for (const auto& b : v)
	    {
	    int	ab_symm { Symm::list[a.item] (Symm::list[b.item]) } ;
	    prod.emplace_back (ab_symm, a.coeff * b.coeff) ;
	    }
	}
    return prod ;
    } ;

SymmSum SymmSum::operator*(SymmTerm b)				// Multiply sum by term
    {
    SymmSum prod ;
    for (const auto& a : *this)
	{
	int ab_symm { Symm::list[a.item] (Symm::list[b.item]) } ;
	prod.emplace_back (ab_symm, a.coeff * b.coeff) ;
	}
    return prod ;
    }

int Symm::operator() (const Symm& a) const noexcept		// Compose tranformations
    {
    Symm c(a) ;
    for (symb x(0) ; x < Nsymb ; ++x)
	{
	c.map[x] = map[a.map[x]] ;
	c.sgn[x] = sgn[a.map[x]] ^ a.sgn[x] ;
	}
    c.sgn[Cbit] = sgn[Cbit] ^ a.sgn[Cbit] ;
    c.name = "" ;
    return trans2indx.at(c) ;
    }

pair<int,Op> Symm::operator()(const Op& a) const		// Transform Op
    {
    Op   b	{ a } ;
    bool neg	{ false } ;

    if (isCodd())
	{
	bool moveE = nonlink(b.front()) && !isferm(b.front()) ;
	if (moveE) rotate(b.begin(), b.begin() + 1, b.end()) ; 
	reverse(b.begin(), b.end()) ;
	}
    for (auto& symb : b)
	{
	neg ^= sgn[symb] ;
	symb = map[symb] ;
	}
    if (!theory.euclid && isferm(b.back()) && isstag(b.back()))
	{
	b.back()  = stag(b.back()) ;
	b.front() = stag(b.front()) ;
	if (b.size() % 2) neg ^= 1 ;
	}
    else if (a.type == OpType::Loop) b.findstart() ;

    return make_pair ((neg ? -1 : 1), b) ;
    }

bool Symm::operator==(const Symm&s) const noexcept		// Compare transformation
    {
    return sgn == s.sgn && std::equal(begin(map), std::end(map), std::begin(s.map)) ;
    }

SymmTerm Symm::known (const string&& name)			// Return named symmetry
    {
    return SymmTerm { list.map.at(name) } ;
    }

string Symm::print() const					// Print transformation
    {
    ostringstream buf ;
    buf << *this ;
    return buf.str() ;
    }

ostream& operator<< (ostream& stream, const Symm& a)		// Print transformation
    {
    int j(0) ;
    stream << a.name << (a.isCodd() ? ":\tC" : ":") ;
    for (symb c(0) ; c < Nsymb ; ++c)
	{
	if (in_thy(c) && (a.map[c] != c || a.sgn[c]))
	    {
	    stream << ((j++%8) ? ", " : "\n\t") << symbname[c] << "->"
		   << (a.sgn[c] ? "-" : "")     << symbname[a.map[c]] ;
	    }
	}
    return stream << "\n" ;
    }

void Symm::symminit()				// Initialize symmetry group data
    {
    uint	indx   (0) ;
    int		dim    { theory.dim } ;
    int		refls  ( ipow (2,dim) ) ;
    symb	link[4]{ '\0', '\1', '\2', '\3' } ;
    char	linkname[4]{ 'x', 'y', 'z', 'w' } ;

    for (int C(0) ; C < 2 ; ++C)			// Loop over C
	{
	vector<int> perm{ 0, 1, 2, 3 } ;

	do  {						// Loop over axis permutations
	    Coord box { theory.box.comp[perm[0]],
			theory.box.comp[perm[1]],
			theory.box.comp[perm[2]],
			theory.box.comp[perm[3]] } ;

	    if (box.point != theory.box.point ) continue ;

	    char	buf[4]{ 'x', 'y', 'z', 'w' } ;
	    string	Pname ;

	    for (int a(0) ; a < dim ; ++a)		// do cyclic decomp of perm
		{
		if (!buf[a] || perm[a] == a) continue ;
		Pname += "P" ;
		Pname += buf[a] ;
		for (int b(perm[a]) ; b != a ; b = perm[b])
		    {
		    Pname += buf[b] ;
		    buf[b] = '\0' ;
		    }
		}

	    for (int R(0) ; R < refls ; ++R)		// Loop over reflections
		{
		Symm	trans ;
		string	CRname ;

		for (symb c(0) ; c < YMend ; ++c)	// Construct symbol mapping
		    {
		    int a {axis(c)} ;
		    trans.map[c] = (c & TandR) | link[perm[a]] ;
		    if (R & (0x1 << perm[a]))
			{
			trans.map[c] ^= Rflag ;
			if (EorElink(c)) trans.sgn.flip(c) ;
			}
		    if (C)
			{
			if (!isE(c) && !isEE(c)) trans.map[c] ^= Rflag ;
			if (EorElink(c)) trans.sgn.flip(c) ;
			}
		    }

		if (C)					// Add conjuation of fermions
		    {
		    for (symb c(YMend) ; c < Nsymb ; ++c)
			{
			if (!theory.euclid)
			    {
			    if (isstag(c) ^ isconj(c)) trans.sgn.flip(c) ;
			    trans.map[c] = stag(refl(c)) ;
			    }
			else
			    {
			    if (isstag(c)) trans.sgn.flip(c) ;
			    trans.map[c] = refl(c) ;
			    }
			}
		    CRname += "C" ;
		    trans.sgn.set(Cbit) ;
		    }
		if (R)					// Compose transfmation name
		    {
		    CRname += "R" ;
		    for (int a(0) ; a < dim ; ++a)
			{
			if (R & (0x1 << a)) CRname += linkname[a] ;
			}
		    }
		trans.indx = indx ;
		trans.name = CRname + Pname ;
		if (trans.name == "") trans.name = "E" ;

		list.push_back(trans) ;
		list.map.insert   (std::make_pair(trans.name, indx)) ;
		trans2indx.insert (std::make_pair(std::move(trans), indx)) ;
		++indx ;
		}
	    } while (next_permutation (perm.begin(), perm.begin() + dim) ) ;
	}
    }
