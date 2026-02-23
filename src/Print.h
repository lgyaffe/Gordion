#ifndef PRINT_H
#define PRINT_H
#include "Gordion.h"

namespace Print
    {
    void print_obs (uint,uint) ;
    void print_obs (uint) ;
    void print_obs (string) ;
    void print_obs () ;
    void print_obs_select (string) ;

    void print_op (uint,uint) ;
    void print_op (uint) ;
    void print_op () ;

    void print_base () ;
    void print_primary () ;

    void print_gen (uint,bool) ;
    void print_gen (bool) ;

    void print_rep (const string&) ;
    void print_rep () ;

    void print_symm (const string&) ;
    void print_symm () ;

    void print_grad (uint,uint) ;
    void print_grad (uint) ;
    void print_grad () ;

    void print_curv (uint,uint,uint) ;
    void print_curv (uint,uint) ;
    void print_curv (uint) ;
    void print_curv () ;

    void print_mode (uint) ;
    void print_mode () ;

    void print_hamiltonian () ;
    void print_freeenergy  () ;
    void print_spectrum    () ;

    void print_geodesic (uint,uint) ;
    void print_geodesic (uint) ;
    void print_geodesic () ;

    void print_lagrange (uint, uint) ;
    void print_lagrange (uint) ;
    void print_lagrange () ;

    void print_blab      () ;
    void print_bcktlist  () ;
    void print_cache     () ;
    void print_couplings () ;
    void print_fermiinit () ;
    void print_geostats	 () ;
    void print_obsstats  () ;
    void print_rkmethods () ;
    void print_state     () ;
    void print_stats     () ;
    void print_symmsets  () ;
    void print_sysindex  () ;
    void print_theory    () ;
    void print_version   () ;
    void print_vevindex  () ;
    } ;

#endif
