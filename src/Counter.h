#ifndef COUNTER_H
#define COUNTER_H
#include "Types.h"
#include "Poly.h"
#include <atomic>
#ifdef DELPHI
#include "Geo.h"
#endif

inline constexpr int cachelinelen = 128 ;	// Cache line size

class alignas (cachelinelen) Counter		// Cache line aligned counter
    {
    std::atomic<ulong> x ;			// Thread-safe

    public:
    Counter& operator++()	 { ++x ; return *this ; }
    Counter& operator=(ulong v)  { x = v ; return *this ; }
    Counter& operator+=(ulong v) { x += v ; return *this ; }
    operator ulong() const	 { return x ; }

    Counter() : x(0) {}

    friend ostream& operator<< (ostream& stream, const Counter& c)
	{ return stream << c.x ; }
    } ;

using CountArr = array<Counter,PSIZ+3> ;	// Array of Counter's

class PolyStats : public CountArr		// Polynomial statistics
    {
    public:
    PolyStats& operator=(const array<ulong,PSIZ+3>& v)
	{
	for (int i(0) ; i < size() ; ++i) (*this)[i] = v[i] ;
	return *this ;
	}
    PolyStats& operator+=(const array<ulong,PSIZ+3>& v)
	{
	for (int i(0) ; i < size() ; ++i) (*this)[i] += v[i] ;
	return *this ;
	}
    void	clear() { std::fill (begin(), end(), 0) ; }
    Counter&	terms() { return (*this)[PSIZ+1] ; }
    Counter&	eqns()  { return (*this)[PSIZ+2] ; }
    ulong	terms() const { return (*this)[PSIZ+1] ; }
    ulong	eqns()  const { return (*this)[PSIZ+2] ; }

    void add (const PolyRec& p)
	{
	array<ulong,PSIZ+3> count {} ;
	auto	recptr	{ p.begin() } ;
	auto	items	{ p.entry().items() } ;
	for (int i(0) ; i < items ; ++i, ++recptr)
	    {
	    const PolyElem*	ptr { &(*recptr) } ;
	    const PolyElem*	end { ptr++->end() } ;
	    doub		val (0.0) ;
	    while (ptr < end)
		{
		PolyTerm t { PolyElem::nextterm (ptr) } ;
		++count[t.order()] ;
		++count[PSIZ+1] ;
		}
	    }
	count[PSIZ+2] += items ;
	*this += count ;
	}

#ifndef DELPHI
#define PW 12
#else
#define PW 9

    void add (const GeoRec& p)
	{
	array<ulong,PSIZ+3> count {} ;
	auto	recptr	{ p.begin() } ;
	auto	items	{ p.entry().items() } ;
	for (int i(0) ; i < items ; ++i, ++recptr)
	    {
	    const GeoElem*	pptr { &(*recptr) } ;
	    const GeoTerm*	ptr { pptr->begin() } ;
	    const GeoTerm*	end { pptr->end() } ;
	    doub		val (0.0) ;
	    while (ptr < end)
		{
		const auto& t { *ptr++ } ;
		++count[t.order()] ;
		++count[PSIZ+1] ;
		}
	    }
	count[PSIZ+2] += items ;
	*this += count ;
	}
#endif

    friend ostream& operator<< (ostream& stream, const PolyStats& s)
	{
	stream << std::setw(PW) << s.eqns() ;
	stream << std::setw(PW) << s.terms() ;
	stream << std::setw(PW-1) << s[0] ;
	for (int i(1) ; i <= PSIZ ; ++i) stream << std::setw(PW) << s[i] ;
	return stream << "\n" ;
	}
    } ;

struct Counters					// Thread-safe counters
    {
    Counter	canons ;	// # Obs canonicalizations
    Counter	symmtrans ;	// # Symmetry transforms
    Counter	commutes ;	// # Commutations
    Counter	stored ;	// # Obs stored
#ifdef DELPHI
    Counter	obs_hit ;	// # Canon::obscache hits
    Counter	obs_miss ;	// # Canon::obscache misses
    Counter	str_hit ;	// # Canon::strcache hits
    Counter	str_miss ;	// # Canon::strcache misses
    Counter	traces ;	// # matrix traces
    Counter	products ;	// # matrix products
    Counter	transforms ;	// # matrix transforms
    Counter	cached ;	// # matrices cached
#else
    Counter	cachehits ;	// # Canon::cache hits
    Counter	cachemiss ;	// # Canon::cache misses
    Counter	assessed ;	// # Obs assessed
    Counter	found ;		// # Obs already known
    Counter	bounded ;	// # Obs xorder bounded
    Counter	discarded ;	// # Obs discarded
    Counter	classified ;	// # Obs classified
    Counter	factored ;	// # Obs factorized
    Counter	reduced ;	// # Obs E/F-reduced
    Counter	approxed ;	// # Obs approximated
    Counter	noapprox ;	// # Obs cannot approximate
    Counter	nointersect ;	// # Obs w/o self-intersection
#endif
    PolyStats	gradstats ;	// # gradient stats
    PolyStats	curvstats ;	// # curvature stats
    PolyStats	lagrstats ;	// # Lagrange brkt stats
    PolyStats	geostats  ;	// # geodesic eqn stats
    } ;

#endif
