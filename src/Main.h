#ifndef MAIN_H
#define MAIN_H
#include "Types.h"

#ifdef DELPHI
#define PRGM delphi
#else
#define PRGM gordion
#endif

#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x)

static const string progname { STRINGIFY(PRGM) } ;
inline const string cmdargs { " [-f <startup_file>] [sys-info or data files]\n" } ;
inline string	    cmdname ;

[[noreturn]] void	quit	  (int) ;
void			sig_catch (int) ;

#endif
