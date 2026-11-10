#include "Blab.h"
#include "Gripe.h"

void Blab::setblab (string file, uint level)	// Set specified blab level
    {
    auto ptr { blabmap.find (file) } ;
    if (ptr != blabmap.end()) blablevel[ptr->second] = level ;
    else gripe ("Unknown file " + file) ;
    }

void Blab::resetblab ()				// (Re)set blab levels
    {
    std::fill (blablevel.begin(), blablevel.end(), 0) ;
    //blablevel[CANON]		= 3 ;
    //blablevel[COMMUTE]	= 3 ;
    //blablevel[GEN]		= 1 ;
    //blablevel[OBS]		= 2 ;
    //blablevel[NUMERICS]	= 1 ;
    }
