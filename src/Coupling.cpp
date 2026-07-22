#include "Coupling.h"
#include "Global.h"

uint Coupling::ncoup ()		// # couplings for stage
    {
    return ncoup (global.stage) ;
    }

uint Coupling::ncoup (int stage)
    {
    return stage ? list.size() : ncoupG ;
    }

string Coupling::values ()	// Make comma-separated coupling values
    {
    std::stringstream buf ;
    string sep {""} ;
    for (const auto& c : Coupling::list)
	{
	if (c.stage > global.stage) continue ;
	buf << sep << c.value ;
	sep = ',' ;
	}
    return string { buf.str() } ;
    }

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
	    if (!e)			continue ;
	    if (e < 0)			stream << "1/" ;
	    				stream << Coupling::list[c].data() ;
	    if (std::abs(e) != 1)	stream << "^" << std::abs(e) ;
					stream << " " ;
	    }
	stream << "( " << term.poly << " ) " ;
	sep = "+ " ;
	}
    if (!vec.front().cpoly.size()) return stream ;
    sep = "" ;
    cout << "\n      = " ;
    for (auto& term : vec)
	{
	stream << sep ;
	for (auto& f : term.coeff)
	    {
	    int  c { f.indx } ;
	    doub e { f.exp  } ;
	    if (!e)			continue ;
	    if (e < 0)			stream << "1/" ;
	    				stream << Coupling::list[c].data() ;
	    if (std::abs(e) != 1)	stream << "^" << std::abs(e) ;
					stream << " " ;
	    }
	stream << "( " << term.cpoly << " ) " ;
	sep = "+ " ;
	}
    return stream ;
    }
