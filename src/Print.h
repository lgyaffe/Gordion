#ifndef PRINT_H
#define PRINT_H
#include "Gordion.h"

namespace Print
    {
    inline bool is_one (doub c) { return std::abs(c - 1) < 1.e-15 ; }
    inline bool is_int (doub c) { return std::abs(c - std::round(c)) < 1.e-15 ; }

    ostream& coeffprt (ostream&, doub) ;

    void print_obs (numb,numb) ;
    void print_obs (numb) ;
    void print_obs () ;
    void print_obs (const string&) ;
    void print_obs_select (const string&) ;

    void print_op (numb,numb) ;
    void print_op (numb) ;
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

    void print_ham_or_free () ;
    void print_spectrum    () ;

    void print_geodesic (numb,uint) ;
    void print_geodesic (numb) ;
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
    void print_stage     () ;
    void print_state     () ;
    void print_stats     () ;
    void print_symmsets  () ;
    void print_sysindex  () ;
    void print_theory    () ;
    void print_version   () ;
    void print_vevindex  () ;
    } ;

#endif
