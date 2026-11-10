#include "Canon.h"
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
	global.baseobs.obsinit	(stage) ;
	Theory::theorydefn	(stage) ;
	OpList::opinit		(stage) ;
	}
    numerics.rk = RKdef() ;
    }
