#include "Global.h"
#include "Numerics.h"

void Global::stageinit (uint i = 0)		// Stage initialization
    {
    if (i && !theory.nf) gripe ("No fermions!!!") ;
    stage = i ? Global::Fermi : Global::Gauge ;
    if (!okfermivev) numerics.initialize (1) ;
    }

string Global::dfltfilename (const string&& ext)  // Make default file name
    {
    const char*	name   { theory.name.data() } ;
    int		approx { approx } ;
    int		obsord { info().maxord } ;
    int		maxgen { info().maxgen } ;

    return format ("{}{}{}{}a{}.{}",
		name, maxgen, fg(), obsord, approx, ext) ;
    }

void Global::mk_bcktlist (uint stage)		// Make Obs bucket list
    {
    auto&	bckt  { info(stage).bckt } ;
    uint	count { ObsList::obs.nobs(stage) } ;
    uint	start { stage ? ObsList::obs.nobsG : 0 } ;
    uint	end   { start + count } ;
    uint	chunk ( 1024 ) ;

    if (count >= MAXBCKT * chunk) chunk = 1 + (count - 1)/MAXBCKT ;

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
    const auto& bckt  { info(stage).bckt } ;

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

void Global::clearpolys (int stage)		// Clear polynomial scripts
    {
    data(stage).grad.clear() ;
    for (auto& mtx  : data(stage).geos) mtx.clear() ;
    for (auto& mtx  : data(stage).lagr) mtx.clear() ;
    for (auto& cube : data(stage).curv) cube.clear() ;
    }
