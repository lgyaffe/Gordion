#include "Hash.h"
#include "Symm.h"

std::size_t Symmhash::operator()(const Symm& s) const noexcept	// Symm hash function
    {
    string buf { s.map.data(), s.map.size() } ;
    std::size_t h1 = std::hash<std::string>{}(buf) ;
    std::size_t h2 = std::hash<std::bitset<Nsymb+1>>{}(s.sgn) ;
    return h1 ^ (h2 << 1) ;
    }

bool Symm_eq::operator()(const Symm& s, const Symm& t) const	// Equality?
    {
    return s.map == t.map && s.sgn == t.sgn ;
    }

std::size_t Strhash::operator()(const Str& s) const		// Str hash function
    {
    return std::hash<string>{}(static_cast<const string&>(s)) ;
    }

