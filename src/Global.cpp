#include "Global.h"
#include "Rep.h"
#include "Canon.h"
#include "Numerics.h"
#include "Gripe.h"

void Global::stageinit (uint i)		// Stage initialization
    {
    if (i == Global::Fermi && !theory.nf) gripe ("No fermions!!!") ;
    stage = i ? Global::Fermi : Global::Gauge ;
    if (!ObsList::base.nobs(stage))	ObsList::base.obsinit () ;
    if (!ObsList::obs.nobs(stage))	ObsList::obs.obsinit  () ;
    if (!info().gens[0].size())		Gen::geninit () ;
    if (!info().Hterms.size())		Theory::theorydefn () ;
    if (i || !numerics.vev.size())	Numerics::numericsinit () ;
    }

string Global::dfltfilename (const string&& ext)  // Make default file name
    {
    const char*	name   { theory.name.data() } ;
    int		approx { approx } ;
    int		obsord { info().maxord } ;
    int		maxgen { info().maxgen } ;

    return format ("{}{}{}{}a{}.{}",
		name, maxgen, global.fg(), obsord, approx, ext) ;
    }

void Global::mk_bcktlist (uint stage)		// Make Obs bucket list
    {
    auto&	info  { global.info(stage) } ;
    auto&	bckt  { info.bckt } ;
    uint	count { ObsList::obs.nobs(stage) } ;
    uint	start { stage ? ObsList::obs.nobsG : 0 } ;
    uint	end   { start + count } ;
    uint	chunk ( 1024 ) ;

    if (count >= MAXBCKT * chunk) chunk = 1 + (count - 1)/MAXBCKT ;
    else while (1 + (count - 1)/chunk > MAXBCKT) chunk *= 2 ;

    bckt.resize (1 + (count - 1)/chunk) ;
    for (uint i(start), k(0) ; i < end ; i += chunk, ++k)
	{
	bckt[k] = uint3 { k, i, std::min(i+chunk, end) - 1 } ;
	}
    }

uint3 Global::bckt_pos (uint i)			// Return stage/bucket/indx
    {
    uint	stage { i >= ObsList::obs.nobsG } ;
    uint	start { stage ? ObsList::obs.nobsG : 0 } ;
    const auto& bckt  { global.info(stage).bckt } ;

    if (bckt.size())
	{
	if (i <= bckt.back()[2])
	    {
	    uint n { bckt.front()[2] - bckt.front()[1] + 1 } ;
	    uint j { (i - start) / n } ;
	    uint k {  i - bckt[j][1]  } ;
	    return uint3 { stage, j, k } ;
	    }
	else fatal ("Invalid observable number!") ;
	}
    else gripe (format ("Make stage {} observables first!", stage)) ;
    }
