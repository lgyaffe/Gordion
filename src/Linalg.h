#ifndef LINALG
#define LINALG
#include "Gordion.h"

#define ARMA_WARN_LEVEL 1
#include <armadillo>

using Uvec = arma::Col<arma::uword> ;
using Rvec = arma::Col<real>  ;
using Dvec = arma::Col<doub>  ;
using Cvec = arma::Col<cmplx> ;
using Dmtx = arma::Mat<doub>  ;
using Cmtx = arma::Mat<cmplx> ;

inline auto opts { arma::solve_opts::refine + arma::solve_opts::equilibrate } ;

#endif
