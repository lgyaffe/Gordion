#include "Canon.h"
#include "Gen.h"
#include "Numerics.h"
#include "Rep.h"
#include "Blab.h"

Global     global ;					// Global information
Numerics   numerics ;					// Numerical data

ObsList	ObsList::obs  {"Canonical",true,true} ;		// Canonical Obs
ObsList	ObsList::base {"Basic"} ;			// Basic defined Obs
ObsList	ObsList::redu {"Reduction"} ;			// Gen reduction Obs

vector<Proj>	Proj::list ;				// Defined Proj's
Index<Symm>	Symm::list ;				// Defined Symm's
Index<Rep>	Rep::list ;				// Defined Rep's
Couplings	Coupling::list ;			// Defined Coupling's
Symmmap		Symm::trans2indx ;			// Symm hash table

void initialize ()
    {
    Blab::resetblab	  () ;
    Theory::theoryinit	  () ;
    Symm::symminit	  () ;
    Rep::repinit	  () ;
    Canon::looptblinit	  () ;
    Canon::spectblinit	  () ;

    for (int stage(0) ; stage < 2 ; ++stage)
	{
	ObsList::base.obsinit	(stage) ;
	Theory::theorydefn	(stage) ;
	Gen::geninit		(stage) ;
	}
    global.stage = Global::Gauge ;
    numerics.rk  = RKdef() ;
    }
