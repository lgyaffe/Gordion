#ifndef BLAB_H
#define BLAB_H
#include "Gordion.h"
#include <map>

class Blab
    {
    public:

    enum BLAB			// Enum distinguishing major .cpp files
	{			// N.B.: want implicit conversion to int
	ASSESS,
	BUILD,
	CANON,
	COMMUTE,
	EVALUATE,
	GEN,
	NUMERICS,
	OBS,
	ODE,
	POLY,
	SAVE,
	SYMB,
	SYMM,
	_BLABNUM_
	} ;
    using Blabmap = std::map<string,BLAB> ;
    using Blabvec = array<uint,_BLABNUM_> ;

    static void	setblab (string,uint) ;			// set level
    static void	resetblab () ;				// reset levels
    static uint level (enum BLAB file)
		    { return blablevel[file] ; }	// return level

    static inline Blabmap blabmap		// file name -> Blab enum
	{
	{"Assess",	ASSESS},
	{"Build",	BUILD},
	{"Canon",	CANON},
	{"Commute",	COMMUTE},
	{"Gen",		GEN},
	{"Numerics",	NUMERICS},
	{"Obs",		OBS},
	{"Ode",		ODE},
	{"Poly",	POLY},
	{"Save",	SAVE},
	{"Symb",	SYMB},
	{"Symm",	SYMM}
	} ;

    private:
    static inline Blabvec blablevel ;		// verbosity levels
    } ;

#endif
