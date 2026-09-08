#include "Canon.h"
#include "Gen.h"
#include "Numerics.h"
#include "Rep.h"
#include "Blab.h"

void initialize ()
    {
    Blab::resetblab	() ;
    Theory::theoryinit	() ;
    Symm::symminit	() ;
    Rep::repinit	() ;
    Canon::looptblinit	() ;
    Canon::spectblinit	() ;

    for (uint stage(0) ; stage < 2 - !theory.nf ; ++stage)
	{
	global.base.obsinit	(stage) ;
	Theory::theorydefn	(stage) ;
	OpList::opinit		(stage) ;
	}
    global.stage = Global::Gauge ;
    numerics.rk  = RKdef() ;
    }
