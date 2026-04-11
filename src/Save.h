#ifndef SAVE_H
#define SAVE_H
#include "Gordion.h"
#include "Coupling.h"
#include "Data.h"
#include "Version.h"
#include <mutex>

namespace Save
    {
    void save_sys	(string="") ;	// Save sys info
    void save_vev	(string="") ;	// Save vev data
    void load_save	(string) ;	// Load save file
    void load_save	(int,string) ;	// Load vev file
    void load_vev	(int) ;		// Load vev set

    void write_header	(fstream&,int=0,int=0) ; // Write file header
    bool read_header	(fstream&,bool) ;	 // Read file header

    void write_sysindex	() ;		// Save SysIndex
    void read_sysindex	() ;		// Load SysIndex

    void write_op	(int) ;		// Save Op's
    void read_op	(int) ;		// Load Op's

    void write_obs	(int) ;		// Save Obs's
    void read_obs	(int) ;		// Load Obs's

    void write_gen	(int) ;		// Save Gen's
    void read_gen	(int) ;		// Load Gen's

    void write_grad	(int) ;		// Save gradient
    void read_grad	(int) ;		// Load gradient

    void write_curv	(int) ;		// Save curvature
    void read_curv	(int) ;		// Load curvature

    void write_lagr	(int) ;		// Save lagrange
    void read_lagr	(int) ;		// Load lagrange

    void write_geos	(int) ;		// Save geodesics
    void read_geos	(int) ;		// Load geodesics

    void write_geo_bckt	(int) ;		// Save geo bucket
    void read_geo_bckt	(int) ;		// Load geo bucket

    void write_coup	() ;		// Save Couplings
    void write_vev	() ;		// Save Obs vev's

    bool read_coup   (int,Couplings*) ;	// Load Couplings
    void read_vev    (int) ;		// Load Obs vev's

    void write_stat	(int)  ;	// Save Statistics
    void read_stat	(int)  ;	// Load Statistics
    void rewrite_stat	() ;		// Save Statistics

    struct FileHdr			// Save file header
	{
	char8	name ;
	Version	version ;
	ushort	ncoup ;
	uint	nvev ;

	ulong vevsize()   const { return nvev  * sizeof (real) ; }
	ulong coupsize()  const { return ncoup * sizeof (Coupling) ; }
	ulong cvsetsize() const { return coupsize() + vevsize() ; }
	bool is_sysfile() const { return !ncoup && !nvev ; }
	bool is_vevfile() const { return  ncoup ||  nvev ; }
	} ;

    inline static FileHdr	filehdr ;	// Save file header
    inline static string	syspath ;	// Sys info file pathname
    inline static std::mutex	savemutex ;	// Sys info file mutex
    } ;

#endif
