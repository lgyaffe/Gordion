#include "Coupling.h"

ostream& operator<< (ostream& stream, const vector<AdjTerm>& vec)	// Print AdjTerm list
    {
    string	plus { " + " } ;
    string	sep ;
    for (auto& term : vec)
	{
	int  c { term.coupindx } ;
	int  e { term.exponent } ;

	stream << sep ;
	if (e)
	    {
	    if (e == -1)	stream << "1/" ;
	    			stream << Coupling::list[c].data() ;
	    if (abs(e) != 1)	stream << "^" << e ;
				stream << " " ;
	    }
	stream << "( " << term.poly << " ) " ;
	//stream << "= ( " << term.cpoly << " ) " ;
	sep = "+ " ;
	}
    return stream ;
    }
