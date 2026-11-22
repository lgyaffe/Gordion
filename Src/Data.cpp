#include "Data.h"
#include "Gen.h"

OpHdr::OpHdr (const Op& op)			// OpHdr constructor
    {
    ulong length { op.size() / sizeof (Element) } ; 
    if (op.size() % sizeof (Element)) ++length ;
    len  = length ;
    type = (int) op.type ;
#ifndef DELPHI
    prim = op.primary ;
    order = op.order ;
    if (order  >= (1 << ordbits)) abort ("Op order exceeds OpHdr capacity") ;
    if (length >= (1 << lenbits)) abort ("Op length exceeds OpHdr capacity") ;
#endif
    }

ObsHdr::ObsHdr (const Obs& obs)			// ObsHdr constructor
    {
    ulong length { obs.size() / sizeof (Element) } ; 
    if (obs.size() % sizeof (Element)) ++length ;
    len  = length ;
    type = (int) obs.type ;
#ifndef DELPHI
    corder = obs.corder ;
    xorder = obs.xorder ;
    if (corder >= (1 << ordbits)) abort ("Obs corder exceeds ObsHdr capacity") ;
    if (xorder >= (1 << ordbits)) abort ("Obs xorder exceeds ObsHdr capacity") ;
    if (length >= (1 << lenbits)) abort ("Obs length exceeds ObsHdr capacity") ;
#endif
    }

GenHdr::GenHdr (const Gen& gen)			// GenHdr constructor
    {
    ulong length { 1 + gen.size() * sizeof (OpTerm) / sizeof (Element) } ;
    len   = length ;
    type  = (int) gen.type ;
    T_odd = gen.T_odd ;
#ifndef DELPHI
    order = gen.order ;
    if (order >= (1 << ordbits))  abort ("Gen order exceeds GenHdr capacity") ;
#else
    group = gen.group ;
    if (group  >= (1 << ordbits))  abort ("Gen group exceeds GenHdr capacity") ;
#endif
    if (length >= (1L << genbits)) abort ("Gen length exceeds GenHdr capacity") ;
    }

