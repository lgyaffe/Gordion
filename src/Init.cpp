#include "Canon.h"
#include "Gen.h"
#include "Numerics.h"
#include "Rep.h"
#include "Blab.h"

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
