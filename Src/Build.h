#ifndef BUILD_H
#define BUILD_H
#include "Op.h"

class Gen ;
class PolyMap ;

namespace Build
    {
    void mk_ham		() ;
    void mk_grad	() ;
    void mk_geos	() ;
    void mk_eqns	(uint) ;
    void mk_curv	(uint) ;
    void mk_lagr	(uint) ;
    void mk_eqns	(string) ;
    void mk_curv	(string) ;
    void mk_lagr	(string) ;

    void mk_obs		(uint) ;
    void mk_loops	() ;
    void mk_Eloops	() ;
    void mk_EEloops	() ;
    void mk_fermions	() ;
    void mk_efermions	() ;

    void do_Loop_bckt	  (const numb3&) ;
    void do_Eloop_bckt	  (const numb3&) ;
    void do_EEloop_bckt   (const numb3&) ;
    void do_Fermion_bckt  (const numb3&) ;
    void do_efermion_bckt (const numb3&) ;
    void do_geo_bckt	  (const numb3&) ;
    void do_geostat_bckt  (const numb3&) ;
    void do_geostats	  () ;

    void check_xorder	(numb, const Gen&, const PolyMap&) ;
    void clear_obs	(uint) ;
    void added_obs	() ;
    } ;

#endif
