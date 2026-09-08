#ifndef STR_H
#define STR_H
#include "Symb.h"

class Str : public string				// Symbol string
    {
    public:
    Str	()		   : string()    {}
    Str	(const symb c)     : string(1,c) {}
    Str	(const symb*    p, const symb*   q) : string(p,q) {}
    Str	(const_iterator p,const_iterator q) : string(p,q) {}
    Str	(const string& s) ;

    string	print (const_iterator, const_iterator) const ;
    string	print () const ;			// Printable version
    int		join (symb) ;				// Append Symb
    int		join (const_iterator, const_iterator) ;	// Append
    int		joinends (Str::iterator, Str::iterator) noexcept ;
    void	excise   (Str::iterator, Str::iterator) noexcept ;	
    bool	isclosed (const_iterator, const_iterator) const noexcept ;
    bool	isclosed () const { return isclosed (cbegin(), cend()) ; }
    Str&	joinends () { resize (joinends (begin(), end())) ; return *this ; }

    static inline bool	dots { false } ;	// print w. symb separators?

    friend ostream& operator<< (ostream&, const Str&) ;
    } ;

struct strhash					// string hash function 
    {
    std::size_t operator()(const string& s) const
	{
	return std::hash<string>{}(s) ;
	}
    using is_transparent = void ;
    } ;

struct str_eq					// string equality function 
    {
    bool operator()(const string& s, const string& t) const
	{
	return std::equal_to<string>{}(s,t) ;
	}
    using is_transparent = void ;
    } ;

#endif
