#include "Coupling.h"

ostream& operator<< (ostream& stream, const vector<AdjTerm>& vec)	// Print AdjTerm list
    {
    string	plus { " + " } ;
    string	sep ;
    for (auto& term : vec)
	{
	stream << sep ;
	for (auto& f : term.coeff)
	    {
	    int  c { f.indx } ;
	    doub e { f.exp  } ;
	    if (!e)		continue ;
	    if (e < 0)		stream << "1/" ;
	    			stream << Coupling::list[c].data() ;
	    if (abs(e) != 1)	stream << "^" << abs(e) ;
				stream << " " ;
	    }
	stream << "( " << term.poly << " ) " ;
	//stream << "= ( " << term.cpoly << " ) " ;
	sep = "+ " ;
	}
    return stream ;
    }
