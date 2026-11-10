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

inline string nicetime (float msecs)
    {
    constexpr sv fmt	{ "{:8.2f} {}" };
    float	 sec	{ 1000 } ;
    float	 min	{ 60 * sec } ;
    float	 hr	{ 60 * min } ;
    float	 day	{ 24 * hr  } ;
    if      (msecs > day)	return format(fmt, msecs / day, "day") ;
    else if (msecs > hr)	return format(fmt, msecs / hr,  "hr")  ;
    else if (msecs > min)	return format(fmt, msecs / min, "min") ;
    else if (msecs > 1)		return format(fmt, msecs / sec, "sec") ;
    else			return format(fmt, msecs,      "msec") ;
    }

inline string nicemem (float rss)
    {
    constexpr sv  fmt	{ "{:8.2f} {}" };
    constexpr int KB	{ 1024 } ;
    constexpr int MB	{ 1024 * KB } ;
    constexpr int GB	{ 1024 * MB } ;
    if      (rss > GB)	return format(fmt, rss / GB, "GB") ;
    else if (rss > MB)	return format(fmt, rss / MB, "MB")  ;
    else if (rss > KB)	return format(fmt, rss / KB, "KB") ;
    else		return format(fmt, rss,      "bytes") ;
    }

#endif
