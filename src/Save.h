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
    string addsubdir	(string,int) ;	// Add subdirectory

    void write_header	(fstream&,uint=0,uint=0) ; // Write file header
    int  read_header	(fstream&) ;	// Read file header

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

    void write_geo_bckt	(int) ;		// Save geo bucket
    void read_geo_bckt	(int) ;		// Load geo bucket

    void write_coup	() ;		// Save Couplings
    void write_vev	() ;		// Save Obs vev's

    bool read_coup   (int,Couplings*) ;	// Load Couplings
    void read_vev    (int) ;		// Load Obs vev's

    void read_stat	() ;		// Load Statistics
    void write_stat	() ;		// Save Statistics
    void rewrite_stat	() ;		// Save Statistics

    ulong vevsize	() ;		// Vev record size
    ulong coupsize	() ;		// Coupling record size
    ulong cvsetsize	() ;		// Coup + Vev size

    static inline struct FileHdr	// Save file header
	{
	char8	name ;
	Version	version ;
	ushort	ncoup ;
	uint	nvev ;
	ulong	hashG ;			// Gauge obs hash
	ulong	hashF ;			// Fermi obs hash

	bool is_sysfile() const { return !ncoup && !nvev ; }
	bool is_vevfile() const { return  ncoup ||  nvev ; }
	} filehdr ;

    inline static std::mutex	savemutex ;	// Sys info file mutex
    } ;

#endif
