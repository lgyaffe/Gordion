#ifndef LIST_H
#define LIST_H
#include "Hash.h"

class Str ;

template <typename T>
class Index : public vector<T>		// vector<T> plus string keyed index
    {
    public:
    hash<string,uint> map ;			// string key -> index map

    numb store (const string& key, const T& t)	// store if new, return indx
	{
	auto [iter, isnew] { map.try_emplace (key, (*this).size()) } ;
	if (isnew) vector<T>::push_back (t) ;
	return iter->second ;
	}
    void clear ()				// Clear index'd list
	{
	map.clear () ;
	vector<T>::clear () ;
	}
    } ;

template <typename T>
class List : protected vector<const T*>		// Indexed list of T*'s
    {
    protected:
    using ListMap = hash<T,numb,Strhash,Str_eq> ;
    ListMap	map ;					// Hash table

    public:
    using const_iterator = typename std::vector<const T*>::const_iterator ;

    const_iterator	cbegin() const noexcept
			    { return vector<const T*>::cbegin() ; }
    const_iterator	cend() const noexcept
			    { return vector<const T*>::cend() ; }

    auto	size() const				// List length
		    { return vector<const T*>::size() ; }
    const T&	operator() (numb indx) const		// Return indexed T
		    { return *(*this)[indx] ; }

    numb	find (const Str& s) const		// Find T, return index
		    {
		    auto p { map.find(s) } ;
		    if (p != map.end()) return p->second ;
		    return MAXNUM ;
		    }
    void	reserve (int len)			// Reserve space
		    {
		    vector<const T*>::reserve (len) ;
		    map.reserve (len) ;
		    }
    } ;

#endif
