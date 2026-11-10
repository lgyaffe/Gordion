#include "Global.h"
#include "Numerics.h"

void Global::stageinit (uint i)			// Stage initialization
    {
    if (i && !theory.nf) gripe ("No fermions!!!") ;
    stage = i ? Global::Fermi : Global::Gauge ;
    }

void Global::close_streams (uint keep)		// Close output streams
    {
    for (uint stage (keep & 1) ; stage < 2 - !theory.nf ; ++stage)
	{
	if (!(keep & 2))  { global.info(stage).sysfile.stream.close() ;
			    global.info(stage).vevfile.stream.close() ; }
	if (!(keep & 4))    global.info(stage).MMAfile.stream.close() ;
	}
    }

void Global::clearpolys (uint stage)		// Clear polynomial scripts
    {
    data(stage).grad.clear() ;
    for (auto& rec : data(stage).geos) rec.clear() ;
    for (auto& rec : data(stage).lagr) rec.clear() ;
    for (auto& rec : data(stage).curv) rec.clear() ;
    }

void Global::mk_bcktlist ()			// Make Obs bucket list
    {
    auto&	bckt	{ info().bckt } ;
    ObsList&	list	{ global.obs } ;
    numb	count	{ list.nobs(stage) } ;
    numb	start	{ stage ? list.nobsG : 0 } ;
    numb	end	( start + count ) ;
    numb	chunk	( 1024 ) ;

    if (count >= MAXBCKT * chunk) chunk = 1 + (count - 1)/MAXBCKT ;

    bckt.resize (count ? 1 + (count - 1)/chunk : 0) ;
    for (numb i(start), k(0) ; i < end ; i += chunk, ++k)
	{
	bckt[k] = numb3 { k, i, std::min(i+chunk, end) - 1 } ;
	}
    }

numb3 Global::bckt_pos (numb i)			// Return stage/bucket/indx
    {
    int		stage { i >= obs.nobsG } ;
    long	start { stage ? obs.nobsG : 0 } ;
    const auto& bckt  { info(stage).bckt } ;

    if (bckt.size())
	{
	if (i <= bckt.back()[2])
	    {
	    numb n { bckt.front()[2] - bckt.front()[1] + 1 } ;
	    numb j ( (i - start) / n ) ;
	    numb k {  i - bckt[j][1]  } ;
	    return numb3 { stage, j, k } ;
	    }
	else fatal ("Invalid observable number!") ;
	}
    else gripe (format("Make stage {} observables first!", stage)) ;
    }

