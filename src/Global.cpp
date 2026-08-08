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

string Global::addsubdir (string dir, int stag)	// Add save sub-directory
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
    for (int stage (keep & 1) ; stage < 2 - !theory.nf ; ++stage)
	{
	if (!(keep & 2))  { global.info(stage).sysfile.stream.close() ;
			    global.info(stage).vevfile.stream.close() ; }
	if (!(keep & 4))    global.info(stage).MMAfile.stream.close() ;
	}
    }

string Global::stageabbrev (int st, const string& ext)
    {
    int		obsord	{ info(st).maxord } ;
    int		genord	{ info(st).maxgen } ;
    string	answer	{ format ("_{}{}{}", genord, fg(st), obsord) } ;
    if (ext != "m")
	{
	ulong	obshash	{ ext == "vev" ? numerics.obshash(st) : info(st).obshash } ;
	string	suffix	{ char('A' + obshash % 26), char('A' + (obshash * 7) % 26) } ;
	answer += suffix ;
	}
    return answer ;
    }

void Global::mk_bcktlist ()			// Make Obs bucket list
    {
    auto&	bckt  { info().bckt } ;
    long	count { info().nobs } ;
    long	start { stage ? info(0).nobs : 0 } ;
    numb	end   ( start + count ) ;
    numb	chunk ( 1024 ) ;

    if (count >= MAXBCKT * chunk) chunk = 1 + (count - 1)/MAXBCKT ;

    bckt.resize (count ? 1 + (count - 1)/chunk : 0) ;
    for (numb i(start), k(0) ; i < end ; i += chunk, ++k)
	{
	bckt[k] = numb3 { k, i, std::min(i+chunk, end) - 1 } ;
	}
    }

numb3 Global::bckt_pos (numb i)			// Return stage/bucket/indx
    {
    int		stage { i >= info(0).nobs } ;
    long	start { stage ? info(0).nobs : 0 } ;
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
    else gripe (format ("Make stage {} observables first!", stage)) ;
    }

void Global::clearpolys (int stage)		// Clear polynomial scripts
    {
    data(stage).grad.clear() ;
    for (auto& rec : data(stage).geos) rec.clear() ;
    for (auto& rec : data(stage).lagr) rec.clear() ;
    for (auto& rec : data(stage).curv) rec.clear() ;
    }
