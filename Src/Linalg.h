#ifndef LINALG
#define LINALG
#include "Types.h"

#define ARMA_WARN_LEVEL 1
#include <armadillo>

using Uvec   = arma::Col<arma::uword> ;
using Rvec   = arma::Col<real>    ;
using Dvec   = arma::Col<doub>    ;
using Cdvec  = arma::Col<dcmplx>  ;
using Dmtx   = arma::Mat<doub>    ;
using Cdmtx  = arma::Mat<dcmplx>  ;

using Cvec   = arma::Col<cmplx>   ;
using Crow   = arma::Row<cmplx>   ;
using Cmtx   = arma::Mat<cmplx>   ;
using Blkmtx = arma::Cube<cmplx>  ;

inline constexpr auto nofill   { arma::fill::none } ;
inline constexpr auto onefill  { arma::fill::ones } ;
inline constexpr auto eyefill  { arma::fill::eye  } ;
inline auto armaopts { arma::solve_opts::refine + arma::solve_opts::equilibrate } ;
#endif
