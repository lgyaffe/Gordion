#ifndef INDEX_H
#define INDEX_H
#include "Gordion.h"

template <typename T>
class Index : public vector<T>		// vector<T> plus index with string keys
    {
    public:
    hash<string,uint>  map ;			// string key -> index map

    uint store (const string& key, const T& t)	// store if new, return indx
	{
	auto [iter, isnew] { map.try_emplace (key, (*this).size()) } ;
	if (isnew) vector<T>::push_back (t) ;
	return iter->second ;
	}
    } ;

#endif
