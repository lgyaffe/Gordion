#ifndef GLOBAL_H
#define GLOBAL_H
#include "Counter.h"
#include "Coupling.h"
#include "Data.h"
#include "Gen.h"
#include "Version.h"
#include <atomic>

#ifdef PARALLEL
#include <tbb/task_arena.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for_each.h>
#define FOR_EACH tbb::parallel_for_each
#define TASK_ARENA(n,cap,code) \
    if (n) { tbb::task_arena arena(n); arena.execute([cap]{code;}); } else code ;
#else
#include <execution>
#define FOR_EACH std::for_each
#define TASK_ARENA(n,cap,code) code
#endif

struct StageInfo				// Stage-specific info
    {
    short			maxgen { 0 } ;	// Max Gen order
    short			maxord { 0 } ;	// Max Obs sc order
    uint			MMAlimit ;	// initial # Obs for MMA file
    ObsSet			MMAlist ;	// addl Obs for MMA file
    OpList			ops ;		// Operators
    array<vector<Gen>,NREP>	gens ;	 	// Generators
    array<ushort,NREP>		neven ;	 	// # T-even generators
    vector<AdjTerm>		Hterms ; 	// Hamiltonian/free energy
    vector<uint3>		bckt ;	 	// Obs bucket list
    Counters			count ;		// Statistics counters
    } ;

static constexpr int MAXBCKT = 256 ;		// Max # Obs buckets

struct SerialData				// Serialised stage data
    {
    DataRec		op   { RecordID::Op }   ; // Operators
    DataRec		obs  { RecordID::Obs }  ; // Observables
    DataRec		gen  { RecordID::Gen }  ; // Generators
    PolyRec		grad { RecordID::Grad } ; // Gradient
    PolyArr<NREP>	curv { RecordID::Curv } ; // Curvature
    PolyArr<NREP>	lagr { RecordID::Lagr } ; // Lagrange bracket
    DataRec		stat { RecordID::Stat } ; // Statistics counters
    PolyArr<MAXBCKT>	geos { RecordID::Geos } ; // Geodesic equations
    } ;

using SysIndex = RecIndxArr<2 * (5 + 2 * NREP + MAXBCKT) + 1> ;
using atombool = std::atomic<bool> ;

class Global					// Global data
    {
    public:
    enum Stage { Gauge = 0, Fermi = 1 } stage ;	// Minimization stage
    SysIndex	sysindex ;			// System data index
    SerialData	stagedata  [2] ;		// Serialized stage data
    StageInfo	stageinfo  [2] ;		// Stage-specific info
    short	repnum     { 0 } ;		// Active irrep number
    short	approx     { 0 } ;		// Approximate Obs's?
    uint	maxthread  { 0 } ;		// Thread limit
    bool	autosave   { false } ;		// Write savefile on bulid
    bool	geoswap    { false } ;		// Swap geo bckts to disk
    bool	massreinit { true } ;		// Auto reinit on new mass?
    bool	symcurv    { true } ;		// Symmetrize curvature?
    bool	fermivev   { false } ;		// Non-base fermi vev's?
    bool	oknegeig   { false } ;		// Negative curvature OK?
    bool	vevappend  { false } ;		// Append to vev data file?
    bool	MMAappend  { false } ;		// Append to MMA file?
    atombool	interrupt  { false } ;		// Interrupt flag
    string	MMAdir     { "./MMA/"  } ;	// MMA result directory
    string	savedir    { "./save/" } ;	// Save file directory
    string	sysfile    ;			// Sys info file
    string	vevfile    ;			// Vev data file
    string	MMAfile    ;			// MMA result file
    fstream	sysstream  ;			// Sys info file stream
    fstream	vevstream  ;			// Vev data file stream
    ofstream	MMAstream  ;			// MMA file stream
    Version	version    ;			// Program version

    auto&	info	(int i)	{ return stageinfo[i] ; }
    auto&	data	(int i)	{ return stagedata[i] ; }
    auto&	info	()	{ return stageinfo[stage] ; }
    auto&	data	()	{ return stagedata[stage] ; }
    auto&	maxgen	()	{ return info().maxgen ; }
    auto&	maxord	()	{ return info().maxord ; }
    auto&	count	()	{ return info().count  ; }
    uint	nobs	()	{ return ObsList::obs.nobs (stage) ; }

    string	dfltfilename (const string&&) ;	// Default file names
    uint3	bckt_pos     (uint) ;		// Obs bucket position
    void	mk_bcktlist  (uint) ;		// Make bucket list
    void	clearpolys   (int) ;		// Clear polys
    void	stageinit    (uint) ;		// Stage initialization
    char	fg (int stage)	const { return stage ? 'f' : 'g' ; }
    char	fg ()		const { return fg (this->stage) ; }

    string sysfilename ()			// Sys info file name
	{
	return sysfile.size() ? sysfile : dfltfilename("sys") ;
	}
    string vevfilename ()			// Vev data file name
	{
	return vevfile.size() ? vevfile : dfltfilename("vev") ;
	}
    string MMAfilename ()			// MMA result file  name
	{
	return MMAfile.size() ? MMAfile : dfltfilename("m") ;
	}
    void streamclose ()				// Close output streams
	{
	sysstream.close() ;
	vevstream.close() ;
	MMAstream.close() ;
	}
    } ;

extern Global global ;


inline DataRec::DataRec (RecordID id)		// DataRec Constructor
    : indexref {global.sysindex.next()}
    { entry().id = id ; }

template <size_t N>
inline PolyArr<N>::PolyArr (RecordID id)	// PolyArr Constructor
    {
    for (auto& p : *this)
	{
	p.indexref = global.sysindex.next()  ;
	p.entry().id = id ;
	}
    }

inline Global global ;				// Global information

#endif
