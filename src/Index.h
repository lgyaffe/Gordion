#ifndef INDEX_H
#define INDEX_H
#include "Gordion.h"

template <typename T>
class Index : public vector<T>		// vector<T> plus index map for
    {					// types T inheriting from string
    public:
    hash<string,uint>  map ;			// name -> index map

    uint store (const string& key, const T& t)	// store if new, return indx
	{
	auto [iter, isnew] { map.try_emplace (key, (*this).size()) } ;
	if (isnew) vector<T>::push_back (t) ;
	return iter->second ;
	}
    uint store (const T& t) { return store (t,t) ; }
    } ;

#endif
