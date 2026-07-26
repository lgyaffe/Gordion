#ifndef TEST_H
#define TEST_H
#include "Gordion.h"

namespace Test
    {
    void irreps() ;					// Test irrep projector orthonormality
    void jacobi() ;					// Test Jacobi identity for Gen triples
    void jacobi(numb) ;					// Test jacobi identities with given Obs
    void jacobi(const string&,const string&,numb) ;	// Test Jacobi identity of given Ops & Obs
    } ;

#endif
