#ifndef SAVE_H
#define SAVE_H
#include "Types.h"
#include "Coupling.h"
#include "Version.h"

namespace Save
    {
    void save_sys	() ;		// Save sys info
    void load_save	(int,string) ;	// Load save file

    string mk_filename	(const string&&) ; // Output file names
    string find_file (const string&, const string&) ;

    int  read_header	(fstream&,const string&); // Read header
    void write_syshead	(fstream&) ; 	// Write sys file header

    void write_sysindex	() ;		// Save SysIndex
    void read_sysindex	() ;		// Load SysIndex

    void write_op	() ;		// Save Op's
    void read_op	() ;		// Load Op's

    void write_obs	() ;		// Save Obs's
    void read_obs	() ;		// Load Obs's

    void write_gen	() ;		// Save Gen's
    void read_gen	() ;		// Load Gen's

    void write_ham	() ;		// Save Hamiltonian
    void read_ham	() ;		// Load Hamiltonian

    void write_grad	() ;		// Save gradient
    void read_grad	() ;		// Load gradient

    void write_curv	() ;		// Save curvature
    void read_curv	() ;		// Load curvature

    void write_lagr	() ;		// Save lagrange
    void read_lagr	() ;		// Load lagrange

    void write_geos	() ;		// Save geodesics
    void read_geos	() ;		// Load geodesics

    void write_coup	 () ;		// Save Couplings
    Couplings* read_coup (int,bool) ;	// Load Couplings

    void read_stat	() ;		// Load Statistics
    void write_stat	() ;		// Save Statistics

    long coupsize	() ;		// Coupling record size
    long datasetsize	() ;		// Total dataset size

    void save_vev	() ;		// Save vev data
    void load_vev	(int) ;		// Load vev set

    void write_vevhead	(fstream&) ;	// Write vev file header
    void write_vev	() ;		// Save Vev's
    void read_vev	(int) ;		// Load Vev's
    long vevsize	() ;		// Vev record size

    string stageabbrev	(uint, const string&) ; // File name info

    void reload_obs	() ;		// Reload Obs
    void rewrite_stat	() ;		// Save Statistics

    void write_geo_bckt	(uint) ;	// Save geo bucket
    void read_geo_bckt	(uint,int) ;	// Load geo bucket

    struct FileHdr			// Save file header
	{
	char8	name ;			// theory name
	Version	version ;		// program version
	ushort	ncoup ;			// # couplings
	union
	    {
	    ulong	nvev ;		// # vev's
	    uint	maxord ;	// max Obs order
	    } ;
	ulong	hashG ;			// Gauge obs hash
	ulong	hashF ;			// Fermion obs hash

	bool is_sysfile() const { return !ncoup ; }
	bool is_vevfile() const { return  ncoup ; }
	} ;

    inline bool	autosave { false } ;	// Auto save flag
    } ;

#endif
