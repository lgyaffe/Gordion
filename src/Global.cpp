#include "Global.h"
#include "Numerics.h"
#include <filesystem>

void Global::stageinit (uint i)			// Stage initialization
    {
    if (i && !theory.nf) gripe ("No fermions!!!") ;
    stage = i ? Global::Fermi : Global::Gauge ;
    }

string Global::mk_filename (const string&& ext)	// Make data file names
    {
    auto	thyname	{ global.stage ? theory.name : theory.parent() } ;
    string	filenam { thyname.data() } ;
    		filenam += stageabbrev(0, ext) ;
    if (stage)	filenam += stageabbrev(1, ext) ;
    if (approx)	filenam += format ("a{}", approx) ;
		filenam += format (".{}", ext) ;
    return filenam ;
    }

string Global::addsubdir (string dir, int stag)	// Return save sub-directory
    {
    dir.append (stag ? theory.name.data() : theory.parent().data() ) ;
    dir += '/' ;
    if (std::filesystem::create_directory (dir))
	cout << "Created subdirectory " << dir << "\n" ;
    return dir ;
    }

string Global::addsubdir (string dir)		// Reteurn save sub-directory
    {
    return addsubdir (dir, stage) ;
    }

void Global::close_streams (int keep)		// Close output streams
    {
    for (int stage (keep & 1) ; stage < 2 ; ++stage)
	{
	if (!(keep & 2))  { global.info(stage).sysstream.close() ;
			    global.info(stage).vevstream.close() ; }
	if (!(keep & 4))    global.info(stage).MMAstream.close() ;
	}
    }

string Global::stageabbrev (int stage, const string& ext)
    {
    int		obsord	{ info(stage).maxord } ;
    int		genord	{ info(stage).maxgen } ;
    string	answer	{ format ("_{}{}{}", genord, fg(stage), obsord) } ;
    if (ext != "m")
	{
	ulong	obshash	{ info(stage).obshash } ;
	string	suffix	{ char('A' + obshash % 26), char('A' + (obshash * 7) % 26) } ;
	answer += suffix ;
	}
    return answer ;
    }

void Global::mk_bcktlist ()			// Make Obs bucket list
    {
    auto&	bckt  { info(stage).bckt } ;
    uint	count { info(stage).nobs } ;
    uint	start { stage ? info(0).nobs : 0 } ;
    uint	end   { start + count } ;
    uint	chunk ( 1024 ) ;

    if (count >= MAXBCKT * chunk) chunk = 1 + (count - 1)/MAXBCKT ;

    bckt.resize (count ? 1 + (count - 1)/chunk : 0) ;
    for (uint i(start), k(0) ; i < end ; i += chunk, ++k)
	{
	bckt[k] = uint3 { k, i, std::min(i+chunk, end) - 1 } ;
	}
    }

uint3 Global::bckt_pos (uint i)			// Return stage/bucket/indx
    {
    uint	stage { i >= info(0).nobs } ;
    uint	start { stage ? info(0).nobs : 0 } ;
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
    for (auto& rec : data(stage).geos) rec.clear() ;
    for (auto& rec : data(stage).lagr) rec.clear() ;
    for (auto& rec : data(stage).curv) rec.clear() ;
    }
