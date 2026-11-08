#ifndef BLAB_H
#define BLAB_H
#include "Types.h"
#include <map>

class Blab
    {
    public:
    enum BLAB			// Enum distinguishing major .cpp files
	{			// N.B.: want implicit conversion to int
	BUILD,
	CANON,
	COMMUTE,
	GEN,
	NUMERICS,
	OBS,
	ODE,
	POLY,
	SAVE,
	STR,
	SYMB,
	SYMM,
#ifdef DELPHI
	LRU,
	MASTER,
	MTX,
#else
	ASSESS,
#endif
	_BLABNUM_
	} ;
    using Blabmap = std::map<string,BLAB> ;
    using Blabvec = array<uint,_BLABNUM_> ;

    static inline Blabmap blabmap		// file name -> Blab enum
	{
	{"Build",	BUILD},
	{"Canon",	CANON},
	{"Commute",	COMMUTE},
	{"Gen",		GEN},
	{"Numerics",	NUMERICS},
	{"Obs",		OBS},
	{"Ode",		ODE},
	{"Poly",	POLY},
	{"Save",	SAVE},
	{"Str",		STR},
	{"Symb",	SYMB},
	{"Symm",	SYMM},
#ifdef DELPHI
	{"LRU",		LRU},
	{"Master",	MASTER},
	{"Mtx",		MTX},
#else
	{"Assess",	ASSESS}
#endif
	} ;

    static void	setblab (string,uint) ;			// set level
    static void	resetblab () ;				// reset levels
    static uint level (enum BLAB file)
		    { return blablevel[file] ; }	// return level

    private:
    static inline Blabvec blablevel ;		// verbosity levels
    } ;

#endif
