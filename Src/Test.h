#ifndef TEST_H
#define TEST_H
#include "Types.h"
#ifdef DELPHI
class Mtx ;
#endif

namespace Test
    {
    void irreps() ;				// Test irrep projector orthonormality
    void jacobi() ;				// Test Jacobi for Gen triples
    void jacobi(numb) ;				// Test Jacobi with given Obs
    void jacobi(const string&,const string&,numb) ; // Test Jacobi for given Ops & Obs
#ifdef DELPHI
    void unitarity() ;				// Test link matrix unitarity
    void unitarity(const string&) ;		// Test holonomy unitarity
    void unitarity(const string&, const Mtx&) ;	// Test holonomy unitarity
#endif
    } ;

#endif
