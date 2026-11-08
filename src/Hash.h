#ifndef HASH_H
#define HASH_H
#include "Types.h"

class Symm ;

template <typename Key, typename... Value>
using hash = unordered_map<Key, Value...> ;	// unordered_map alias

template <uint N>
struct Indxhash					// index array hash function
    {
    std::size_t operator() (const array<numb,N>& t) const noexcept
	{
	std::hash<numb> hasher ;
	size_t answer { 0 } ;
	for (auto x : t)
	    {
	    answer ^= hasher(x) + 0x9e3779b9 + (answer << 6) + (answer >> 2) ;
	    }
	return answer;
	}
    using is_transparent = void ;
    } ;

struct Strhash					// string hash function 
    {
    std::size_t operator()(const string& s) const
	{
	return std::hash<string>{}(s) ;
	}
    using is_transparent = void ;
    } ;

struct Str_eq					// string equality function 
    {
    bool operator()(const string& s, const string& t) const
	{
	return std::equal_to<string>{}(s,t) ;
	}
    using is_transparent = void ;
    } ;

struct Symmhash					// Symm hash function 
    {
    std::size_t operator()(const Symm&) const noexcept ;
    } ;

struct Symm_eq					// Symm equality test
    {
    bool operator()(const Symm&, const Symm&) const ;
    } ;

#endif
