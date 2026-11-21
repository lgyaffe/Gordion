#include "Data.h"
#include "Gen.h"


OpHdr::OpHdr (const Op& op)			// OpHdr constructor
    {
    type = (int) op.type ;
    len  = op.size() / sizeof (Element) ; 
    if (op.size() % sizeof (Element)) ++len ;
#ifndef DELPHI
    prim = op.primary ;
    order = op.order ;
    if (order >= (1 << ordbits)) abort ("Op order exceeds OpHdr capacity") ;
    if (len   >= (1 << lenbits)) abort ("Op length exceeds OpHdr capacity") ;
#endif
    }

ObsHdr::ObsHdr (const Obs& obs)			// ObsHdr constructor
    {
    type = (int) obs.type ;
    len  = obs.size() / sizeof (Element) ; 
    if (obs.size() % sizeof (Element)) ++len ;
#ifndef DELPHI
    corder = obs.corder ;
    xorder = obs.xorder ;
    if (corder >= (1 << ordbits)) abort ("Obs corder exceeds ObsHdr capacity") ;
    if (xorder >= (1 << ordbits)) abort ("Obs xorder exceeds ObsHdr capacity") ;
    if (len    >= (1 << lenbits)) abort ("Obs length exceeds ObsHdr capacity") ;
#endif
    }

GenHdr::GenHdr (const Gen& gen, int repnum)	// GenHdr constructor
    {
    len   = 1 + gen.size() * sizeof (OpTerm) / sizeof (Element) ;
    type  = (int) gen.type ;
    T_odd = gen.T_odd ;
    rep   = repnum ;
#ifndef DELPHI
    order = gen.order ;
    if (order >= (1 << ordbits)) abort ("Gen order exceeds GenHdr capacity") ;
#else
    group = gen.group ;
    if (group >= (1 << ordbits)) abort ("Gen group exceeds GenHdr capacity") ;
#endif
    if (len   >= (1 << lenbits)) abort ("Gen length exceeds GenHdr capacity") ;
    }

