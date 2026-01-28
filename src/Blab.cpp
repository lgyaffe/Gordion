#include "Blab.h"
#include "Gripe.h"

void Blab::setblab (string file, uint level)	// Set specified blab level
    {
    auto ptr { blabmap.find (file) } ;
    if (ptr != blabmap.end()) blablevel[ptr->second] = level ;
    else gripe (format ("Unknown file {}",file)) ;
    }

void Blab::resetblab ()				// (Re)set blab levels
    {
    std::fill (blablevel.begin(), blablevel.end(), 0) ;
    //blablevel[BLAB::ASSESS]	= 4 ;
    //blablevel[BLAB::CANON]	= 3 ;
    //blablevel[BLAB::COMMUTE]	= 3 ;
    //blablevel[BLAB::GEN]	= 1 ;
    //blablevel[BLAB::OBS]	= 2 ;
    //blablevel[BLAB::NUMERICS]	= 1 ;
    }
