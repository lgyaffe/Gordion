#ifndef BUILD_H
#define BUILD_H
#include "Op.h"
#include <mutex>

class Gen ;
class PolyMap ;

namespace Build
    {
    void mk_obs		(int) ;
    void mk_ham		() ;
    void mk_grad	() ;
    void mk_curv	(uint) ;
    void mk_lagr	(uint) ;
    void mk_curv	(string) ;
    void mk_lagr	(string) ;
    void mk_geos	() ;
    void mk_loops	() ;
    void mk_Eloops	() ;
    void mk_EEloops	() ;
    void mk_fermions	() ;
    void mk_Efermions	() ;
    void do_geostats	() ;

    void do_geo_bckt	  (const numb3&) ;
    void do_Loop_bckt	  (const numb3&) ;
    void do_Eloop_bckt	  (const numb3&) ;
    void do_EEloop_bckt   (const numb3&) ;
    void do_Fermion_bckt  (const numb3&) ;
    void do_Efermion_bckt (const numb3&) ;
    void do_geostat_bckt  (const numb3&) ;
    void check_xorder	  (numb, const Gen&, const PolyMap&) ;
    void clear_obs	  (int) ;
    void close_streams	  () ;

    inline static int		cord ;			// Current creation order
    inline static Obsset	newobs ;		// Newly generated Obs
    inline static vector<Op>	primepOps ;		// primary conj mom Ops
    inline static std::mutex	obsmutex ;		// newobs insertion lock
    constexpr int		single_thread = 1024 ;	// multi-threading threshold
    } ;

#endif
