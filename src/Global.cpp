#include "Global.h"
#include "Rep.h"
#include "Canon.h"
#include "Numerics.h"
#include "Gripe.h"

void Global::stageinit (uint stage)		// Stage initialization
    {
    if (stage == Global::Fermi && !theory.nf) gripe ("No fermions!!!") ;
    global.stage = stage ? Global::Fermi : Global::Gauge ;
    Global::clearbuild     (true) ;
    ObsList::base.obsinit  () ;
    Gen::geninit	   () ;
    Op::setprimary	   () ;
    Theory::theorydefn	   () ;
    Numerics::numericsinit () ;
    }

void Global::clearbuild (bool zapobs)		// Clear prior build data
    {
    for (uint i (global.stage) ; i < 2 ; ++i)
	{
	global.data(i).grad.clear() ;
	for (auto& mtx  : global.data(i).geos) mtx.clear() ;
	for (auto& mtx  : global.data(i).lagr) mtx.clear() ;
	for (auto& cube : global.data(i).curv) cube.clear() ;
	}
    if (zapobs) clearobs () ;
    }

void Global::clearobs ()			// Clear observables
    {
    if (global.stage == 0)
	{
	Canon::cache.clear() ;
	global.info(0).nobs = 1 ;
	global.info(1).nobs = 0 ;
	global.info(0).maxord = 2 ;
	global.info(1).maxord = 2 ;
	}
    else
	{
	Canon::cache.purge (global.nobsG()) ;
	global.info(1).nobs = 0 ;
	global.info(1).maxord = 2 ;
	}
    ObsList::obs.purge (global.nobsG()) ;
    ObsList::obs.obsinit () ;
    }

void Global::mk_bcktlist (uint stage)		// Make Obs bucket list
    {
    auto&	info  { global.info(stage) } ;
    auto&	bckt  { info.bckt } ;
    uint	count { info.nobs } ;
    uint	start { stage ? global.nobsG() : 0 } ;
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
    uint	stage { i >= global.nobsG() } ;
    uint	start { stage ? global.nobsG() : 0 } ;
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
    else gripe ("Make observables first!") ;
    }

string Global::dfltfilename (const string&& ext)  // Make default file name
    {
    const char*	name   { theory.name.data() } ;
    int		approx { global.approx } ;
    int		obsord { global.info().maxord } ;
    int		maxgen { global.info().maxgen } ;
    int		minlim { numerics.minlim } ;
    int		minord { minlim ? std::min (minlim,maxgen) : maxgen } ;
    int		genord { ext == "sys" ? maxgen : minord } ;

    return format ("{}{}{}{}a{}.{}",
		name, genord, global.fg(), obsord, approx, ext) ;
    }
