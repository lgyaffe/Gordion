#ifndef STR_H
#define STR_H
#include "Symb.h"
#ifdef DELPHI
class Symm ;
#endif

class Str : public string				// Symbol string
    {
    public:
    Str	()		   : string()    {}
    Str	(const symb c)     : string(1,c) {}
    Str	(const symb*    p,const symb*    q) : string(p,q) {}
    Str	(const_iterator p,const_iterator q) : string(p,q) {}
    Str	(const string& s) ;

    string	print (const_iterator, const_iterator) const ;
    string	print    () const ;			// Printable version
    int		join     (symb) ;			// Append Symb
    int		join     (const_iterator, const_iterator) ; // Append

    int		joinends (Str::iterator, Str::iterator) noexcept ;
    void	excise   (Str::iterator, Str::iterator) noexcept ;	
    bool	isclosed (const_iterator, const_iterator) const noexcept ;
    bool	isclosed () const { return isclosed (cbegin(), cend()) ; }
    Str&	joinends () { resize (joinends (begin(), end())) ; return *this ; }

#ifdef DELPHI
    bool	alllinks () const ;			// No specials?
    bool	nongauge () const ;			// Is non-gauge??
    void	validate () const ;			// Validate
    int		canon    () ;				// Canonicalize
    int		trans    (const Symm&) noexcept ;	// Symmetry tranform
#endif

    static bool compare (const Str& a, const Str& b) noexcept { return a < b ; }
    static inline bool	dots { false } ;	// print w. symb separators?

    friend ostream& operator<< (ostream&, const Str&) ;
    } ;

#endif
