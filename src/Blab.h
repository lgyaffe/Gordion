#ifndef BLAB_H
#define BLAB_H
#include "Gordion.h"
#include <map>

enum BLAB			// Enum distinguishing major .cpp files
    {				// N.B.: want implicit conversion to int
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

class Blab
    {
    public:
    using Blabmap = std::map<string,BLAB> ;
    using Blabvec = array<uint,_BLABNUM_> ;

    static void	setblab(string,uint) ;		// set verbosity level
    static void	resetblab() ;			// reset verbosity levels
    static inline Blabvec blablevel ;		// verbosity levels
    static inline Blabmap blabmap		// file name -> Blab enum
	{
	{"Assess",	BLAB::ASSESS},
	{"Build",	BLAB::BUILD},
	{"Canon",	BLAB::CANON},
	{"Commute",	BLAB::COMMUTE},
	{"Gen",		BLAB::GEN},
	{"Numerics",	BLAB::NUMERICS},
	{"Obs",		BLAB::OBS},
	{"Ode",		BLAB::ODE},
	{"Poly",	BLAB::POLY},
	{"Save",	BLAB::SAVE},
	{"Symb",	BLAB::SYMB},
	{"Symm",	BLAB::SYMM}
	} ;
    } ;

#endif
