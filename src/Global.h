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

static constexpr int NENTRY = 6 + 2 * NREP + MAXBCKT + 1 ;

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

struct SaveFile					// Save file info
    {
    string	path ;
    fstream	stream ;
    bool	append ;
    } ;

struct MMAFile					// MMA output file info
    {
    string	path ;
    ofstream	stream ;
    bool	append ;
    ObsSubset	obs ;
    } ;

struct StageInfo				// Stage-specific info
    {
    using Genvec = vector<Gen> ;

    short		maxgen   { 0 } ;	// Max Gen order
    short		maxord   { 0 } ;	// Max Obs sc order
    bool		validvev { false } ;	// Vev list valid?
    OpList		ops ;			// Operators
    array<Genvec,NREP>	gens ;	 		// Generators
    array<ushort,NREP>	neven ;	 		// # T-even generators
    vector<AdjTerm>	Hterms ; 		// Hamiltonian/free energy
    vector<numb3>	bckt ;	 		// Obs bucket list
    SaveFile		sysfile ;		// Sys-info file
    SaveFile		vevfile ;		// Vev-data file
    MMAFile		MMAfile ;		// MMA results file
    SysIndex		sysindex ;		// System data index
    Counters		count ;			// Statistics counters
    } ;

class ObsInfo : public ObsList			// Canonical Obs info
    {
    public:
    using ObsList::ObsList ;

    static inline numb			nobsG ;		// shadows ObsList::nobsG
    static inline numb			nobsF ;		// shadows ObsList::nobsF
    static inline ulong			hashG ;		// shadows ObsList::hashG
    static inline ulong			hashF ;		// shadows ObsList::hashF
    static inline std::vector<numb2>	fermiinit ;	// Fermion -> Loop map
    static inline bool			swapped {false};// Swapped to disk?
    static inline thread_local bool	freeze {true} ;	// Freeze master list?
    static inline thread_local Obsset	inbox ;		// Obs awaiting insertion
    } ;

class Global					// Global data
    {
    public:
    using atombool = std::atomic<bool> ;
    using Stage = enum { Gauge = 0, Fermi = 1 } ;	
    
    StageInfo	stageinfo [2] ;			// Stage-specific info
    SerialData	stagedata [2] { stageinfo[0].sysindex,
				stageinfo[1].sysindex } ;

    ObsList	base {"Basic" } ;		// Basic Obs
    ObsInfo	obs  {"Canonical",true,true} ;	// Canonicalized Obs

    Stage	stage ;				// Minimization stage
    short	repnum     { 0 } ;		// Active irrep number
    short	approx     { 0 } ;		// Approximate Obs's?
    uint	maxthread  { 0 } ;		// Thread limit
    bool	autosave   { false } ;		// Write savefile on bulid
    bool	geoswap    { false } ;		// Swap geo bckts to disk
    bool	obsswap    { false } ;		// Swap obs list to disk
    atombool	interrupt  { false } ;		// Interrupt flag
    string	savedir    { "./save/" } ;	// Save file directory
    string	MMAdir	   { "./MMA/"  } ;	// MMA result directory
    Version	version    ;			// Program version

    auto&	info	(int i)	{ return stageinfo[i] ; }
    auto&	data	(int i)	{ return stagedata[i] ; }
    auto&	info	 ()	{ return stageinfo[stage] ; }
    auto&	data	 ()	{ return stagedata[stage] ; }
    auto&	maxgen	 ()	{ return info().maxgen ; }
    auto&	maxord	 ()	{ return info().maxord ; }
    auto&	count	 ()	{ return info().count  ; }
    numb	nobs	 ()	{ return stage ? obs.nobsF : obs.nobsG ; }

    string	stageabbrev (int, const string&) ; // File name info
    string	mk_filename (const string&&) ;	// Output file names
    string	addsubdir     (string, int) ;	// Add theory subdir
    string	addsubdir     (string) ;	// Add theory subdir
    void	mk_bcktlist   ()    ;		// Make bucket list
    void	close_streams (int) ;		// Close output streams
    void	clearpolys    (int) ;		// Clear polys
    numb3	bckt_pos      (numb) ;		// Obs bucket position
    void	stageinit     (uint) ;		// Stage initialization
    char	fg (int stage)	const { return stage ? 'f' : 'g' ; }
    char	fg ()		const { return fg (this->stage) ; }
    } ;

inline Global global ;					// Global information

inline DataRec::DataRec (SysIndex& indx, RecordID id)	// DataRec constructor
    : indexref { indx.next() }
    { entry().id = id ; }

template <size_t N>					// PolyArr constructor
template <size_t... Is> constexpr
inline PolyArr<N>::PolyArr (std::index_sequence<Is...>, SysIndex& indx, RecordID id)
    : array<PolyRec,N> { PolyRec ( (static_cast<void>(Is), indx), id )... } 
    {}

#endif
