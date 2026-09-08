#ifndef SAVE_H
#define SAVE_H
#include "Gordion.h"
#include "Coupling.h"
#include "Data.h"
#include "Version.h"
#include <mutex>

namespace Save
    {
    void save_sys	() ;		// Save sys info
    void save_vev	() ;		// Save vev data
    void load_save	(string) ;	// Load save file
    void load_save	(int,string) ;	// Load vev file
    void load_vev	(int) ;		// Load vev set
    void reload_obs	() ;		// Reload Obs
    string addsubdir	(string,int) ;	// Add subdirectory

    void write_syshead	(fstream&) ; 	// Write sys file header
    void write_vevhead	(fstream&) ;	// Write vev file header
    int  read_header	(fstream&,const string&) ; // Read file header

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

    void write_geo_bckt	(uint) ;	// Save geo bucket
    void read_geo_bckt	(uint,uint) ;	// Load geo bucket

    void write_coup	 () ;		// Save Couplings
    Couplings* read_coup (int,bool) ;	// Load Couplings

    void write_vev	() ;		// Save Vev's
    void read_vev	(int) ;		// Load Vev's

    void read_stat	() ;		// Load Statistics
    void write_stat	() ;		// Save Statistics
    void rewrite_stat	() ;		// Save Statistics

    long vevsize	() ;		// Vev record size
    long coupsize	() ;		// Coupling record size
    long datasetsize	() ;		// Coup + Vev record size

    static inline struct FileHdr	// Save file header
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
	} filehdr ;

    inline static std::mutex	savemutex ;	// Sys info file mutex
    } ;

#endif
