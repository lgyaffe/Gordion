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

static constexpr int MAXBCKT = 256 ;		// Max # Obs buckets

struct SerialData				// Serialised stage data
    {
    SerialData (SysIndex&) ;		// Constructor

    DataRec		op    ;		// Operators
    DataRec		obs   ;		// Observables
    DataRec		gen   ;		// Generators
    PolyRec		ham   ;		// Hamiltonian
    PolyRec		grad  ;		// Gradient
    PolyArr<NREP>	curv  ;		// Curvature
    PolyArr<NREP>	lagr  ;		// Lagrange bracket
    DataRec		stat  ;		// Statistics counters
    PolyArr<MAXBCKT>	geos  ;		// Geodesic equations
    } ;

static constexpr int NENTRY = 6 + 2 * NREP + MAXBCKT ;

class SysIndex : public RecIndxArr<NENTRY> {} ;	// Sys-info file index

inline SerialData::SerialData (SysIndex& indx)	// SerialData constructor
    :
    op   { indx, RecordID::Op   },
    obs  { indx, RecordID::Obs  },
    gen  { indx, RecordID::Gen  },
    ham  { indx, RecordID::Ham  },
    grad { indx, RecordID::Grad },
    curv { indx, RecordID::Curv },
    lagr { indx, RecordID::Lagr },
    stat { indx, RecordID::Stat },
    geos { indx, RecordID::Geos }
    {}

struct StageInfo				// Stage-specific info
    {
    using Genvec = vector<Gen> ;

    short		maxgen  { 0 } ;		// Max Gen order
    short		maxord  { 0 } ;		// Max Obs sc order
    uint		nobs    { 0 } ;		// # Obs
    ulong		obshash { 0 } ;		// Obs list hash
    OpList		ops ;			// Operators
    array<Genvec,NREP>	gens ;	 		// Generators
    array<ushort,NREP>	neven ;	 		// # T-even generators
    vector<AdjTerm>	Hterms ; 		// Hamiltonian/free energy
    vector<uint3>	bckt ;	 		// Obs bucket list
    string		savedir { "./save/" } ;	// Save file directory
    string		syspath ;		// Sys info path
    fstream		sysstream ; 		// Sys info stream
    string		vevpath ;		// Vev data path
    fstream		vevstream ; 		// Vev data stream
    string		MMAdir  { "./MMA/"  } ;	// MMA result directory
    string		MMApath ;		// MMA results path
    ObsSubset		MMAobs ;		// MMA Obs list
    ofstream		MMAstream ; 		// MMA results stream
    SysIndex		sysindex ;		// System data index
    Counters		count ;			// Statistics counters
    } ;

class Global					// Global data
    {
    public:
    using atombool = std::atomic<bool> ;
    using Stage = enum { Gauge = 0, Fermi = 1 } ;
    
    Stage	stage ;
    StageInfo	stageinfo [2] ;
    SerialData	stagedata [2] { stageinfo[0].sysindex,
				stageinfo[1].sysindex } ;

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
    Version	version    ;			// Program version

    auto&	info	(int i)	{ return stageinfo[i] ; }
    auto&	data	(int i)	{ return stagedata[i] ; }
    auto&	info	()	{ return stageinfo[stage] ; }
    auto&	data	()	{ return stagedata[stage] ; }
    auto&	maxgen	()	{ return info().maxgen ; }
    auto&	maxord	()	{ return info().maxord ; }
    auto&	count	()	{ return info().count  ; }
    uint	nobs	()	{ return info().nobs   ; }

    string	stageabbrev (int, const string&) ; // File name info
    string	mk_filename (const string&&) ;	// Default file names
    void	mk_bcktlist () ;		// Make bucket list
    uint3	bckt_pos    (uint) ;		// Obs bucket position
    void	clearpolys  (int) ;		// Clear polys
    void	stageinit   (uint) ;		// Stage initialization
    char	fg (int stage)	const { return stage ? 'f' : 'g' ; }
    char	fg ()		const { return fg (this->stage) ; }
    } ;

extern Global global ;					// Global information

inline DataRec::DataRec (SysIndex& indx, RecordID id)	// DataRec constructor
    : indexref { indx.next() }
    { entry().id = id ; }

template <size_t N>					// PolyArr constructor
template <size_t... Is> constexpr
inline PolyArr<N>::PolyArr (std::index_sequence<Is...>, SysIndex& indx, RecordID id)
    : std::array<PolyRec,N> { PolyRec ( (static_cast<void>(Is), indx), id )... } 
    {}

#endif
