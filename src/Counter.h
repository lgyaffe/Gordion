#ifndef COUNTER_H
#define COUNTER_H
#include "Theory.h"
#include <atomic>

static constexpr int cachelinelen = 128 ;	// Cache line size

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

template <size_t N>
class CountArr : public array<Counter,N>	// Array of Counter's
    {
    public:
    CountArr& operator=(ulong v)
	{
	for (int i(0) ; i < N ; ++i) (*this)[i] = v ;
	return *this ;
	}
    CountArr& operator=(array<ulong,N> v)
	{
	for (int i(0) ; i < N ; ++i) (*this)[i] = v[i] ;
	return *this ;
	}
    CountArr& operator+=(array<ulong,N> v)
	{
	for (int i(0) ; i < N ; ++i) (*this)[i] += v[i] ;
	return *this ;
	}
    } ;

struct Counters					// Thread-safe counters
    {
    Counter		canons ;	// # Obs canonicalizations
    Counter		cachehits ;	// # Canon::cache hits
    Counter		cachemiss ;	// # Canon::cache misses
    Counter		symmtrans ;	// # Symmetry transforms
    Counter		commutes ;	// # Commutations
    Counter		assessed ;	// # Obs assessed
    Counter		found ;		// # Obs already known
    Counter		bounded ;	// # Obs xorder bounded
    Counter		discarded ;	// # Obs discarded
    Counter		classified ;	// # Obs classified
    Counter		factored ;	// # Obs factorized
    Counter		reduced ;	// # Obs E/F-reduced
    Counter		stored ;	// # Obs stored
    Counter		approxed ;	// # Obs approximated
    Counter		noapprox ;	// # Obs cannot approximate
    Counter		nointersect ;	// # Obs w/o self-intersection
    Counter		ngeos ;		// # geodesic equations
    Counter		geoterms ;	// # geodesic equation terms
    CountArr<PSIZ+1>	geotermord ;	// # monomial terms by order

    void cleargeostats ()		// Clear geodesics counters
	{
	ngeos = 0 ;
	geoterms = 0 ;
	geotermord = 0 ;
	}
    } ;

#endif
