#ifndef BUILD_H
#define BUILD_H
#include "Gordion.h"

class Gen ;
class PolyMap ;

namespace Build
    {
    void mk_obs(int) ;
    void mk_loops() ;
    void mk_Eloops() ;
    void mk_EEloops() ;
    void mk_fermions() ;
    void mk_Efermions() ;
    void mk_geos() ;
    void mk_grad() ;
    void mk_curv(int) ;
    void mk_lagr(int) ;
    void do_geostats() ;
    void do_geo_bckt	  (const uint3&) ;
    void do_Loop_bckt	  (const uint3&) ;
    void do_Eloop_bckt	  (const uint3&) ;
    void do_EEloop_bckt   (const uint3&) ;
    void do_Fermion_bckt  (const uint3&) ;
    void do_Efermion_bckt (const uint3&) ;
    void do_geostat_bckt  (const uint3&) ;
    void check_xorder (uint, const Gen&, const PolyMap&) ;
    } ;

#endif
