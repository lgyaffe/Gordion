#ifndef TERM_H
#define TERM_H
#include "Print.h"

template <typename C,typename T> class Term 	// T object w. C coefficient
    {
    public:
    C	coeff = 0.0 ;
    T	item ;

    Term () : item() {}
    Term (const T t, C d = 1.0) : item(t), coeff(d) {}

    friend ostream& operator<< (ostream& stream, const Term<C,T>& x)
	{
	Print::coeffprt (stream, x.coeff) ;
	return stream << x.item ;
	}
    } ;

#ifdef DELPHI
class IndxSet : public vector<numb2>		// index range list
    {
    numb prevsize = 128 ;			// last compress size
    public:
    bool do_compress() const			// trigger compression?
	    { return size() > 2 * prevsize ; }

    void compress()				// sort & compress
	    {
	    if ((prevsize = size()))
		{
		std::sort (begin(), end()) ;
		auto ptr  { begin() } ;
		auto next { ptr + 1 } ;
		for (; next < end() ; ++next)
		    {
		    if (ptr->second + 1 >= next->first)	// merge entries
			{
			if (ptr->second < next->second)
			    ptr->second = next->second ;
			}
		    else *++ptr = *next ;
		    }
		resize (ptr + 1 - begin()) ;
		}
	    }

    friend ostream& operator<< (ostream& stream, const IndxSet& set)
	{
	string sep { "" } ;
	for (const auto& [b,e]: set)
	    {
	    if (b==e)	stream << sep << b ;
	    else	stream << sep << b << "-" << e ;
	    sep = "," ;
	    }
	return stream ;
	}
    } ;
#endif

#endif
